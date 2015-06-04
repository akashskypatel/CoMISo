//=============================================================================
//
//  CLASS DOCloudSolver - IMPLEMENTATION
//
//=============================================================================

// TODO: this uses Cbc for the MPS file export; consider implementing our own 
// MPS/LP export to remove the dependency.

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_DOCLOUD_AVAILABLE

//=============================================================================
#include "DOCloudSolver.hh"

#include <Base/Debug/DebTime.hh>
#include <Base/Debug/DebUtils.hh>
#include <Base/Utils/OutcomeUtils.hh>

#include <curl/curl.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <thread>         
#include <chrono>         

#include <stdio.h>
#include <stdlib.h>
#include <io.h>

DEB_module("DOCloudSolver")

//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION ==========================================================

namespace {

namespace cURLpp { // some classes to wrap around the libcurl C data

struct Session
{
  Session() { curl_global_init(CURL_GLOBAL_DEFAULT); }
  ~Session() { curl_global_cleanup(); }
};

// TODO: This inheritance could be restrictive in the future for now is OK
class Request
{
public:
  Request() : hnd_(curl_easy_init()), http_hdr_(nullptr) {}
  
  virtual ~Request() 
  { 
    curl_easy_cleanup(hnd_); 
    if (http_hdr_ != nullptr)
      curl_slist_free_all(http_hdr_);
  }

  //operator CURL*() { return hnd_; } 
  //operator const CURL*() const { return hnd_; } 

  bool valid() const { return hnd_ != nullptr; }

  void set_url(const char* _url)
  {
    curl_easy_setopt(hnd_, CURLOPT_URL, _url);
  }

  void add_http_header(const char* _hdr)
  {
    http_hdr_ = curl_slist_append(http_hdr_, _hdr);
  }

  // TODO: REMOVE this security hole!!
  void disable_ssl_verification()
  {
    /*
     * If you want to connect to a site who isn't using a certificate that is
     * signed by one of the certs in the CA bundle you have, you can skip the
     * verification of the server's certificate. This makes the connection
     * A LOT LESS SECURE.
     *
     * If you have a CA cert for the server stored someplace else than in the
     * default bundle, then the CURLOPT_CAPATH option might come handy for
     * you.
     */ 
    curl_easy_setopt(hnd_, CURLOPT_SSL_VERIFYPEER, 0L);
 
    /*
     * If the site you're connecting to uses a different host name that what
     * they have mentioned in their server certificate's commonName (or
     * subjectAltName) fields, libcurl will refuse to connect. You can skip
     * this check, but this will make the connection less secure.
     */ 
    curl_easy_setopt(hnd_, CURLOPT_SSL_VERIFYHOST, 0L);
  }

  void perform()
  {
    DEB_enter_func;

    prepare();

    // set the header we send to the server
    curl_easy_setopt(hnd_, CURLOPT_HTTPHEADER, http_hdr_);
    // set the write function to handle incoming data from the server
    curl_easy_setopt(hnd_, CURLOPT_WRITEFUNCTION, write_func);
    // set the string to store the incoming header data
    curl_easy_setopt(hnd_, CURLOPT_HEADERDATA, reinterpret_cast<void*>(&hdr_));
    // set the string to store the incoming main body (data)
    curl_easy_setopt(hnd_, CURLOPT_WRITEDATA, reinterpret_cast<void*>(&bdy_));
    // do the transmission
    CURLcode res;
    int try_nmbr = 0;
    do
    {
      // CURLE_SSL_CONNECT_ERROR is something that we are seeing a lot with 
      // the DOcloud service, a single retry is usually sufficient to recover.
      DEB_line_if(try_nmbr > 0, 2, "curl_easy_perform() retry #" << try_nmbr);
      res = curl_easy_perform(hnd_);
      if (res != CURLE_OK)
      {
        DEB_warning(1, "curl_easy_perform() failed with code: " << res << 
          ", message: " << curl_easy_strerror(res));
      }
    } while (res == CURLE_SSL_CONNECT_ERROR && try_nmbr++ < 10);

    THROW_OUTCOME_if(res != CURLE_OK, TODO);

    DEB_line(6, "Received Header: " << hdr_);
    DEB_line(6, "Received Body: " << bdy_);

    finalize();
  }

  const std::string& header() const { return hdr_; }
  const std::string& body() const { return bdy_; }

protected:
  static size_t write_func(const char* _ptr, const size_t _size, 
    const size_t _nmemb, void* _str)
  {
    // TODO: not sure how much exception-safe (e.g. out of memory this is!)
    size_t n_add = _size * _nmemb;
    if (n_add == 0)
      return 0;

    auto& str = *reinterpret_cast<std::string*>(_str);
    str.append(_ptr, n_add);
    return n_add;
  }

  // Apparently read callback function is required with Windows libcurl DLL!?
  //static size_t read_func(void* _ptr, const size_t _size, 
  //  const size_t _nmemb, FILE* _file)
  //{
  //  return fread(_ptr, _size, _nmemb, _file);
  //}


protected:
  CURL* hnd_;

protected:
  //virtual functions to control the perform() behavior
  virtual void prepare() {}
  virtual void finalize() {}

private:
  curl_slist* http_hdr_;
  std::string hdr_;
  std::string bdy_;
};

class Post : public Request
{
public: 
  Post(const char* _post) : post_(_post) {}
  Post(const std::string& _post) : post_(_post) {}

protected:
    virtual void prepare()
    {
      // set the post fields
      curl_easy_setopt(hnd_, CURLOPT_POSTFIELDS, post_.data());
      curl_easy_setopt(hnd_, CURLOPT_POSTFIELDSIZE, post_.size());
    }

private:
  std::string post_;
};

class Upload : public Request
{
public:
  Upload(const char* _filename) : filename_(_filename), file_(nullptr) {}
  Upload(const std::string& _filename) : filename_(_filename), file_(nullptr) {}
  virtual ~Upload()
  {
    if (file_ != nullptr)
      std::fclose(file_);
  }

protected:
  virtual void prepare();
  virtual void finalize();

private:
  std::string filename_;
  FILE* file_;
};

void Upload::prepare()
{
  /* tell it to "upload" to the URL */ 
  curl_easy_setopt(hnd_, CURLOPT_UPLOAD, 1L);

  /* set where to read from (on Windows you need to use READFUNCTION too) */ 
  file_ = std::fopen(filename_.data(), "rb");
  curl_easy_setopt(hnd_, CURLOPT_READDATA, file_);

  /* and give the size of the upload (optional) */ 
  const auto filelen = _filelength(fileno(file_));
  curl_easy_setopt(hnd_, CURLOPT_INFILESIZE_LARGE, (curl_off_t)filelen);

  /* we want to use our own read function */
  //curl_easy_setopt(hnd_, CURLOPT_READFUNCTION, read_func);

  /* enable verbose for easier tracing */ 
  //curl_easy_setopt(hnd_, CURLOPT_VERBOSE, 1L);
}

void Upload::finalize()
{
  DEB_enter_func;

  double rate, time;

  /* now extract transfer info */ 
  curl_easy_getinfo(hnd_, CURLINFO_SPEED_UPLOAD, &rate);
  curl_easy_getinfo(hnd_, CURLINFO_TOTAL_TIME, &time);

  DEB_double_format("%.2f");
  DEB_line(2, 
    "Upload speed: " << rate / 1024. << "Kbps; Time: " << time << "s.");
}

class Get : public Request
{
protected:
    virtual void prepare() { curl_easy_setopt(hnd_, CURLOPT_HTTPGET, 1L); }
};

class Delete : public Request
{
protected:
    virtual void prepare() 
    { 
      curl_easy_setopt(hnd_, CURLOPT_CUSTOMREQUEST, "DELETE"); 
    }
};

}//cURLpp 

namespace DOcloud {

static char* root_url__ = 
  "https://api-oaas.docloud.ibmcloud.com/job_manager/rest/v1/jobs";
static std::string api_key__ = 
  "X-IBM-Client-Id: api_0821c92f-0f2b-4ea5-be24-ecc9cd7695dd";
static char* app_type__ = "Content-Type: application/json";

class HeaderTokens
{
public:
  HeaderTokens(const std::string& _hdr)
  {
    // TODO: Performance can be improved by indexing, strtok_r(), etc ...
    //  ... but probably not worth the effort
    std::istringstream strm(_hdr);
    typedef std::istream_iterator<std::string> Iter;
    std::copy(Iter(strm), Iter(), std::back_inserter(tkns_));  
  }

  const std::string& operator[](const size_t _idx) const
  {
    return tkns_[_idx];
  }

  size_t number() const { return tkns_.size(); }

  // Find a token equal to the label and return its value (next token)
  bool find_value(const std::string& _lbl, std::string& _val) const
  {
    auto it = std::find(tkns_.begin(), tkns_.end(), _lbl);
    if (it == tkns_.end() || ++it == tkns_.end())
      return false;

    _val = *it;
    return true;
  }

  typedef std::vector<std::string>::const_iterator const_iterator;

  const_iterator begin() const { return tkns_.begin();}
  const_iterator end() const { return tkns_.end();}

private:
  std::vector<std::string> tkns_;
};

class JsonTokens
{
public:
  JsonTokens() {}
  JsonTokens(const std::string& _bdy) { set(_bdy); }

  void set(const std::string& _bdy)
  {
    ptree_.clear();
    if (_bdy.empty())
      return;
    std::istringstream strm(_bdy);
    boost::property_tree::json_parser::read_json(strm, ptree_);
  }

  //size_t number() const { return tkns_.size(); }

  // Find a token equal to the label and return its value
  bool find_value(const std::string& _lbl, std::string& _val) const
  {
    auto it = ptree_.find(_lbl);
    if (it == ptree_.not_found())
      return false;

    _val = it->second.get_value<std::string>();
    return true;
  }

  typedef boost::property_tree::ptree PTree;

  const PTree& ptree() const { return ptree_; }

private:
  PTree ptree_;
};

Debug::Stream& operator<<(Debug::Stream& _ds, const JsonTokens::PTree& _ptree)
{
  std::stringstream os;
  boost::property_tree::json_parser::write_json(os, _ptree);
  _ds << os.str();
  return _ds;
}

Debug::Stream& operator<<(Debug::Stream& _ds, const JsonTokens& _json_tkns)
{
  return _ds << _json_tkns.ptree();
}

void throw_http_error(const int _err_code, const std::string& _bdy)
{
  DEB_enter_func;

  std::string err_msg;
  JsonTokens bdy_tkns(_bdy);
  bdy_tkns.find_value("message", err_msg);
  DEB_warning(1, "HTTP Status Code: " << _err_code << "; Message: " << err_msg);

  switch (_err_code)
  {
  case 400 : THROW_OUTCOME(TODO); // Invalid job creation data / status
  case 403 : THROW_OUTCOME(TODO); // Subscription limit exceeded
  case 404 : THROW_OUTCOME(TODO); // Requested job could not be found
  default : THROW_OUTCOME(TODO); // Unrecognized HTTP status code
  }
}

void check_http_error(
  const cURLpp::Request& _rqst, 
  const HeaderTokens& _hdr_tkns,
  const int code_ok = 201
  )
{
  const std::string http_lbl = "HTTP/1.1";
  const int code_cntn = 100; // continue code, ignore
  
  for (auto it = _hdr_tkns.begin(), it_end = _hdr_tkns.end(); it != it_end; 
    ++it)
  {
    if (*it != http_lbl) // search for the http label token
      continue; 
    THROW_OUTCOME_if(++it == it_end, TODO); // missing http code
    const auto code = atoi(it->data());

    if (code == code_ok) 
      return; // success code found, exit here
    else if (code == code_cntn)
      continue; // continue code found, continue
    else // another code found, throw an error
      throw_http_error(code, _rqst.body());
  }
  THROW_OUTCOME(TODO); // final http code not found
}


class Job : public cURLpp::Session
{
public:
  Job(const char* _filename) : filename_(_filename) {}
  Job(const std::string& _filename) : filename_(_filename) {}
  ~Job();

  void setup()
  {
    make();
    upload();
    start();
  }

  void wait();
  void sync_status();
  void sync_log();
  bool active() const; // requires synchronized status
  bool stalled() const 
  { 
    // exit quick if we have a solution, or wait 5 min if we don't have one
    return (sol_nmbr_ > 0 && stld_sec_nmbr_ >= 15) || 
      (sol_nmbr_ == 0 && stld_sec_nmbr_ >= 300);
  }

  void abort();
  void solution(std::vector<double>& _x) const;

private:
  const std::string filename_;
  std::string url_;
  JsonTokens stts_;
  // these variables are initialized in start() 
  int log_seq_idx_; // the log sequence number, used to get DOcloud log entires
  int sol_nmbr_; // number of solutions found so far, according to the log
  int sol_sec_nmbr_; // number of seconds at the last new solution 
  int stld_sec_nmbr_; // number of seconds since the last new solution 

private:
  void make();
  void upload();
  void start();
};

Job::~Job()
{
  DEB_enter_func;

  if (url_.empty()) // not setup
    return;

  cURLpp::Delete del;
  if (!del.valid()) 
  {
    DEB_error("Failed to construct a delete request");
    return; // no point in throwing an exception here
  }

  del.set_url(url_.data());
  del.add_http_header(api_key__.c_str());
  del.perform();

  // no point in checking the return value either, we can't do much if the
  // delete request has failed
}

void Job::make()
{
  DEB_enter_func;

  cURLpp::Post post(std::string(
    "{\"attachments\" : [{\"name\" :\"" + filename_ + "\"}]}"));
  THROW_OUTCOME_if(!post.valid(), TODO); //Failed to initialize the request

  post.set_url(root_url__);
  post.add_http_header(api_key__.c_str());
  post.add_http_header(app_type__);
  post.perform();

  HeaderTokens hdr_tkns(post.header());
  check_http_error(post, hdr_tkns);
  // TODO: DOcloud header is successful but no location value
  THROW_OUTCOME_if(!hdr_tkns.find_value("Location:", url_), TODO); 
}

void Job::upload()
{
  cURLpp::Upload upload(filename_);
  THROW_OUTCOME_if(!upload.valid(), TODO); //Failed to initialize the request

  auto url = url_ + "/attachments/" + filename_ + "/blob";
  upload.set_url(url.data());
  upload.add_http_header(api_key__.c_str());
  upload.perform();
  HeaderTokens hdr_tkns(upload.header());
  check_http_error(upload, hdr_tkns, 204);
}

void Job::start()
{
  cURLpp::Post post("");
  THROW_OUTCOME_if(!post.valid(), TODO); //Failed to initialize the request

  auto url = url_ + "/execute";
  post.set_url(url.data());
  post.add_http_header(api_key__.c_str());
  post.add_http_header(app_type__);
  post.perform();
  HeaderTokens hdr_tkns(post.header());
  check_http_error(post, hdr_tkns, 204);

  log_seq_idx_ = sol_nmbr_ = sol_sec_nmbr_ = stld_sec_nmbr_ = 0;
}

void Job::sync_status()
{
  DEB_enter_func;

  cURLpp::Get get;
  THROW_OUTCOME_if(!get.valid(), TODO); //Failed to initialize the request
  get.set_url(url_.data());
  get.add_http_header(api_key__.c_str());
  get.perform();
  HeaderTokens hdr_tkns(get.header());
  check_http_error(get, hdr_tkns, 200);

  stts_.set(get.body());
/*
  // The code below attempted to analyse the status data to find out the 
  // progress of the solver. This is an undocumented use and does not seem to 
  // work so far. Achieved here for potential use in the future.
  
  DEB_line(2, stts_);

  const auto& details = stts_.ptree().get_child_optional("details");
  if (!details)
    return;

  DEB_line(2, details.get());

  const auto& prg_gap = 
    details.get().get_child("PROGRESS_GAP").get_value<std::string>();
  std::string mip_gap;
  const auto mip_gap_it = details.get().find("cplex.mipabsgap");
  if (mip_gap_it != details.get().not_found())
    mip_gap = mip_gap_it->second.get_value<std::string>();

  DEB_line(2, "Status, MIP gap: " << mip_gap << "; Progress gap: " << prg_gap);
*/
}

void Job::sync_log()
{
  DEB_enter_func;

  cURLpp::Get get;
  THROW_OUTCOME_if(!get.valid(), TODO); //Failed to initialize the request
  const std::string url = url_ + "/log/items?start=" + 
    std::to_string(log_seq_idx_) + "&continuous=true";
  get.set_url(url.data());
  get.add_http_header(api_key__.c_str());
  get.perform();
  HeaderTokens hdr_tkns(get.header());
  check_http_error(get, hdr_tkns, 200);

  JsonTokens log(get.body());

  // iterate the log items, deb_out messages and analyze for solutions #
  for (const auto& log_item : log.ptree())
  {
    DEB_line_if(log_seq_idx_ == 0, 2, "**** DOcloud log ****");
    const auto& records = log_item.second.get_child("records");
    for (const auto& record : records)
    {// the message ends with \n
      const std::string msg = record.second.get_child("message").
        get_value<std::string>();
      DEB_out(2, record.second.get_child("level").get_value<std::string>() << 
        ": " << msg);

      const int time_str_len = 15;
      const char time_str[time_str_len + 1] = "Elapsed time = ";
      const auto time_str_idx = msg.find(time_str);
      if (time_str_idx == std::string::npos)
        continue;

      const int sec_nmbr = atoi(msg.data() + time_str_idx + time_str_len);
      //DEB_line(1, "# seconds elapsed : " << sec_nmbr);

      const int sol_str_len = 12;
      const char sol_str[sol_str_len + 1] = "solutions = ";
      const auto sol_str_idx = msg.find(sol_str);
      if (sol_str_idx == std::string::npos)
        continue;

      const int sol_nmbr = atoi(msg.data() + sol_str_idx + sol_str_len);
      //DEB_line(1, "# solutions found so far: " << sol_nmbr);
      if (sol_nmbr > sol_nmbr_) // new solution(s) found
      {// update the number of solutions and the time of the last solution found
        sol_nmbr_ = sol_nmbr;
        sol_sec_nmbr_ = sec_nmbr;
      }
      stld_sec_nmbr_ = sec_nmbr - sol_sec_nmbr_;
    }
    log_seq_idx_ = log_item.second.get_child("seqid").get_value<int>() + 1; 
  }
}

bool Job::active() const
{
  std::string exct_stts;
  stts_.find_value("executionStatus", exct_stts);
  
  // assume the job is not active if the status is not recognized
  return exct_stts == "CREATED" || exct_stts == "NOT_STARTED" || 
    exct_stts == "RUNNING" || exct_stts == "INTERRUPTING";

  /*
  Backup of old code converting execution status strings to enum value
  enum StatusType { ST_CREATED, ST_NOT_STARTED, ST_RUNNING, ST_INTERRUPTING, 
    ST_INTERRUPTED, ST_FAILED, ST_PROCESSED, ST_UNKNOWN };

  const int n_stts = (int)ST_UNKNOWN;
  const char stts_tbl[n_stts][16] = { "CREATED", "NOT_STARTED", "RUNNING", 
    "INTERRUPTING", "INTERRUPTED", "FAILED", "PROCESSED" };

  for (int i = 0; i < n_stts; ++i)
  {
    if (stts == stts_tbl[i])
      return (StatusType)i;
  }
  return ST_UNKNOWN;
  while (stts == ST_CREATED || stts == ST_NOT_STARTED || stts == ST_RUNNING ||
      stts == ST_INTERRUPTING);
  */
}

void Job::abort()
{
  std::string exct_stts;
  stts_.find_value("executionStatus", exct_stts);
  if (exct_stts != "RUNNING")
    return; // already aborted or aborting

  cURLpp::Delete del;
  THROW_OUTCOME_if(!del.valid(), TODO); //Failed to initialize the request
  const std::string url = url_ + "/execute";
  del.set_url(url.data());
  del.add_http_header(api_key__.c_str());
  del.perform();

  HeaderTokens hdr_tkns(del.header());
  check_http_error(del, hdr_tkns, 204);
}

void Job::wait()
{
  do 
  {
    std::this_thread::sleep_for(std::chrono::seconds(1));
    sync_status();
    sync_log();
    if (stalled())
      abort();
  } while (active());
}

void Job::solution(std::vector<double>& _x) const
{
  DEB_enter_func;

  // check the solution status (assume it's synchronized already)

  // What are the possible values for this??
  std::string slv_stts;
  stts_.find_value("solveStatus", slv_stts);
  DEB_line(2, "solveStatus=" << slv_stts);
  
  cURLpp::Get get;
  THROW_OUTCOME_if(!get.valid(), TODO); //Failed to initialize the request

  auto url = url_ + "/attachments/solution.json/blob";
  get.set_url(url.data());
  get.add_http_header(api_key__.c_str());
  get.perform();

  HeaderTokens hdr_tkns(get.header());
  check_http_error(get, hdr_tkns, 200);

  JsonTokens bdy_tkns(get.body());
  DEB_line(7, bdy_tkns);

  const auto& vrbls = bdy_tkns.ptree().get_child("CPLEXSolution.variables");
  const auto n_vrbls = vrbls.size();
  THROW_OUTCOME_if(n_vrbls != _x.size(), TODO); // Solution variables number does not match
  
  size_t i = 0;
  for (const auto& v : vrbls)
  {
    // TODO: this way of conversion is rather hacky
    const std::string name = 
      v.second.get_child("name").get_value<std::string>(); // this is x#IDX
    const int idx = atoi(name.data() + 1);
    THROW_OUTCOME_if(idx < 0 || idx > n_vrbls, TODO);  // Invalid index
    _x[idx] = v.second.get_child("value").get_value<double>();

    DEB_out(7, "#" << idx << "=" << 
      v.second.get_child("value").get_value<std::string>() << "; ");
  }

  DEB_line(3, "X=" << _x);
}


} // namespace DOcloud

std::string lp_file_name()
{
  // TODO: This is not MT-safe, it's not even process-safe!
  static int n_mps = 0;
  std::string filename("DOcloud_problem_");
  filename += std::to_string(n_mps++);
  return filename + ".lp";
}

#define P(X) ((X).data())
#define XVAR(IDX) "x" << IDX

class WriteExpression
{
  // lp format allows a max line length = 560.
  // For simplicity we put end line when the number of written characters on the
  // same line exceeds LINE_TRESHOLD_LEN.
  enum { LINE_TRESHOLD_LEN = 100 };

public:
  WriteExpression(std::ofstream& _out_str) : out_str_(_out_str)
  {
    start();
  }

  void start()
  {
    f_size_ = out_str_.tellp();
    at_start_ = true;
  }
  
  // Writes a monomial.
  void add_monomial(const double _coeff, const size_t _i_var)
  {
    if (_coeff == 0)
      return;
    if (_coeff == 1)
    {
      if (!at_start_)
        out_str_ << " + ";
    }
    else if (_coeff == -1)
      out_str_ << " - ";
    else
    {
      if (!at_start_)
      {
        if (_coeff > 0)
          out_str_ << " + ";
        else
          out_str_ << ' ';
      }
      out_str_ << _coeff << ' ';
    }
    out_str_ << XVAR(_i_var);

    at_start_ = false;

    const auto new_f_size = out_str_.tellp();
    if (new_f_size - f_size_ > LINE_TRESHOLD_LEN)
    {
      out_str_ << std::endl;
      f_size_ = new_f_size;
    }
  }

private:
  std::ofstream& out_str_;
  std::fstream::pos_type f_size_;
  bool at_start_;
};

// Create a lp file for the given constraints and object function.
// Here is the lp format specifications:
// http://www-01.ibm.com/support/knowledgecenter/SSSA5P_12.6.1/ilog.odms.cplex.help/CPLEX/FileFormats/topics/LP.html
std::string create_lp_file(
  NProblemInterface* _problem,
  const std::vector<NConstraintInterface*>& _constraints,
  const std::vector<PairIndexVtype>& _discrete_constraints,
  const std::vector<double>& _x
  )
{
  const int n_cols = _problem->n_unknowns(); // Unknowns #

  const std::string f_name(lp_file_name());
  std::ofstream lp_file(f_name);
  THROW_OUTCOME_if(!lp_file.is_open(), TODO);

  // Set the ofstream options.
  lp_file << std::setprecision(std::numeric_limits<double>::digits10 + 2);

  lp_file << "\\Problem name: " << std::endl << std::endl;
  lp_file << "Minimize" << std::endl;

  // Writes objective function.
  lp_file << "obj: ";

  std::vector<double> objective(n_cols);
  _problem->eval_gradient(P(_x), P(objective));
  WriteExpression wrte_expr(lp_file);
  for (size_t i = 0; i < objective.size(); ++i)
    wrte_expr.add_monomial(objective[i], i);

  // Writes constraints.
  lp_file << std::endl << "Subject To" << std::endl;
  for (const auto& cstr : _constraints)
  {
    NConstraintInterface::SVectorNC gc;
    cstr->eval_gradient(P(_x), gc);

    wrte_expr.start();
    for (NConstraintInterface::SVectorNC::InnerIterator v_it(gc); v_it; ++v_it)
    {
      auto coeff = v_it.value();
      wrte_expr.add_monomial(coeff, v_it.index());
    }
    switch (cstr->constraint_type())
    {
    case NConstraintInterface::NC_EQUAL:
      lp_file << " = ";
      break;
    case NConstraintInterface::NC_GREATER_EQUAL:
      lp_file << " >= ";
      break;
    case NConstraintInterface::NC_LESS_EQUAL:
      lp_file << " <= ";
      break;
    default:
      THROW_OUTCOME(TODO); // Can not express current constraint.
    }
    lp_file << -cstr->eval_constraint(P(_x)) << std::endl;
  }

  // Writes the variables.
  lp_file << "Bounds" << std::endl;
  for (size_t i = 0; i < n_cols; ++i)
    lp_file << XVAR(i) << " Free" << std::endl;

  // Integer and binary variables.
  std::vector<unsigned int> int_var, bin_var;
  for (const auto& dc : _discrete_constraints)
  {
    if (dc.second == Integer)
      int_var.push_back(dc.first);
    else if (dc.second == Binary)
      bin_var.push_back(dc.first);
  }
  auto write_var_set = [&lp_file](const std::vector<unsigned int>& _vars,
    const char* _type)
  {
    if (_vars.empty())
      return;
    // Writes integer variables.
    lp_file << _type << std::endl;
    auto var_it = _vars.begin();
    lp_file << XVAR(*var_it);
    size_t n_wrt_var = 1;
    while (++var_it != _vars.end())
    {
      if (n_wrt_var++ % 16) // 16 variables per line. Lines length must be < 560.
        lp_file << ' ';
      else
        lp_file << std::endl;
      lp_file << XVAR(*var_it);
    }
    lp_file << std::endl;

  };
  // Writes integer variables.
  write_var_set(int_var, "Integers");

  // Writes Binary variables.
  write_var_set(bin_var, "Binary");

  lp_file << "End";

  return f_name;
}

#undef XVAR 

}// namespace 


void DOCloudSolver::set_api_key(const char* _api_key)
{
  DOcloud::api_key__ = std::string("X-IBM-Client-Id: ") + _api_key;
}

void DOCloudSolver::solve(
  NProblemInterface*                        _problem,
  const std::vector<NConstraintInterface*>& _constraints,
  const std::vector<PairIndexVtype>&        _discrete_constraints,
  const double                              _time_limit
)
{
  DEB_enter_func;
  DEB_warning_if(!_problem->constant_hessian(), 1,
    "DOCloudSolver received a problem with non-constant hessian!");
  DEB_warning_if(!_problem->constant_gradient(), 1,
    "DOCloudSolver received a problem with non-constant gradient!");

  std::vector<double> x(_problem->n_unknowns(), 0.0); // solution
  std::string filename = create_lp_file(_problem, _constraints, 
    _discrete_constraints, x);

  DOcloud::Job job(filename);
  job.setup();
  job.wait();
  job.solution(x);
  
  _problem->store_result(P(x));
}

#undef P

//=============================================================================
} // namespace COMISO
//=============================================================================

#endif // COMISO_DOCLOUD_AVAILABLE
//=============================================================================

