//=============================================================================
//
//  namespace COMISO::DOcloud IMPLEMENTATION
//
//=============================================================================

#include "DOcloudJob.hh"
#if COMISO_DOCLOUD_AVAILABLE

#include "DOCloudConfig.hh"
#include <Base/Utils/OutcomeUtils.hh>
#include <Base/Debug/DebUtils.hh>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>
#include <thread>         

DEB_module("DOcloud")
  
namespace COMISO {
namespace DOcloud { 

//////////////////////////////////////////////////////////////////////////
// Config
static const char* app_type__ = "Content-Type: application/json";

Config::Config()
  : root_url_("https://api-oaas.docloud.ibmcloud.com/job_manager/rest/v1/jobs"),
    api_key_("X-IBM-Client-Id: api_0821c92f-0f2b-4ea5-be24-ecc9cd7695dd"),
    infs_time_(300), fsbl_time_(15), 
    cache_loc_("\\\\camfs1\\General_access\\Martin_Marinov\\ReForm\\Cache\\")
{
  const char* env_cache_dir = getenv("ReFormCacheDir");
  if (env_cache_dir != nullptr && env_cache_dir[0] != 0)
  {
    cache_loc_ = env_cache_dir;
    if (cache_loc_.back() != '\\')
      cache_loc_ += '\\'; // Eventually add '\' to the directory string.
  }
}

Config& Config::object()
{
  // TODO: implement MT-lock
  static Config config;
  return config;
}

void Config::set_root_url(const char* const _root_url)
{
  root_url_ = _root_url;
}

void Config::set_api_key(const char* _api_key)
{
  api_key_ = std::string("X-IBM-Client-Id: ") + _api_key;
}

void Config::set_infeasible_timeout(const int _infs_time)
{
  infs_time_ = _infs_time;
}

void Config::set_feasible_timeout(const int _fsbl_time)
{
  fsbl_time_ = _fsbl_time;
}

void Config::set_cache_location(const char* const _cache_loc)
{
  cache_loc_ = _cache_loc;
}

//////////////////////////////////////////////////////////////////////////
// Config
const Config& Config::query() { return object(); } 
Config& Config::modify() { return object(); } 

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

  const_iterator begin() const { return tkns_.begin(); }
  const_iterator end() const { return tkns_.end(); }

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

class HttpStatus
{
public:
  HttpStatus(const cURLpp::Request& _rqst)
    : rqst_(_rqst), hdr_tkns_(_rqst.header()), code_(0)
  {
    const std::string http_lbl = "HTTP/1.1";
    const int code_cntn = 100; // continue code, ignore

    for (auto it = hdr_tkns_.begin(), it_end = hdr_tkns_.end(); it != it_end; ++it)
    {
      if (*it != http_lbl) // search for the http label token
        continue; 
      THROW_OUTCOME_if(++it == it_end, TODO); // missing http code
      code_ = atoi(it->data());
      if (code_ != code_cntn)
        return;
    }
    THROW_OUTCOME(TODO); // final http code not found
  }

  void check(int _code_ok) const
  {
    if (code_ != _code_ok)  // another code found, throw an error
      throw_http_error(code_, rqst_.body());
  }

  const int& code() const { return code_; }

  const HeaderTokens& header_tokens() const { return hdr_tkns_; }

private:
  const cURLpp::Request& rqst_;
  HeaderTokens hdr_tkns_;
  int code_;
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
  del.add_http_header(Config::query().api_key());
  del.perform();

  delete stts_;
  // no point in checking the return value either, we can't do much if the
  // delete request has failed
}

void Job::make()
{
  DEB_enter_func;
  std::string post_loc = 
    "{\"attachments\" : [{\"name\" :\"" + std::string(filename_) + "\"}]}";
  cURLpp::Post post(post_loc);
  THROW_OUTCOME_if(!post.valid(), TODO); //Failed to initialize the request

  post.set_url(Config::query().root_url());
  post.add_http_header(Config::query().api_key());
  post.add_http_header(app_type__);
  post.perform();

  HttpStatus http_stat(post);
  http_stat.check(201);
  // TODO: DOcloud header is successful but no location value
  THROW_OUTCOME_if(!http_stat.header_tokens().find_value("Location:", url_),
    TODO); 

  if (stts_ == nullptr) 
    stts_ = new JsonTokens;
}

void Job::upload(cURLpp::Upload& _upld)
{
  THROW_OUTCOME_if(!_upld.valid(), TODO); //Failed to initialize the request

  auto url = url_ + "/attachments/" + filename_ + "/blob";
  _upld.set_url(url.data());
  _upld.add_http_header(Config::query().api_key());
  _upld.perform();
  HttpStatus http_stat(_upld);
  http_stat.check(204);
}

void Job::upload()
{
  if (file_buf_.empty())
  {// file is not buffered into memory
    cURLpp::UploadFile upld(filename_);
    upload(upld);
  }
  else
  {
    cURLpp::UploadData upld(file_buf_);
    upload(upld);
  }
}

void Job::start()
{
  cURLpp::Post post("");
  THROW_OUTCOME_if(!post.valid(), TODO); //Failed to initialize the request

  auto url = url_ + "/execute";
  post.set_url(url.data());
  post.add_http_header(Config::query().api_key());
  post.add_http_header(app_type__);
  post.perform();
  HttpStatus http_stat(post);
  http_stat.check(204);

  log_seq_idx_ = sol_nmbr_ = sol_sec_nmbr_ = stld_sec_nmbr_ = 0;
  stop_wtch_.start();
}

void Job::sync_status()
{
  DEB_enter_func;

  cURLpp::Get get;
  THROW_OUTCOME_if(!get.valid(), TODO); //Failed to initialize the request
  get.set_url(url_.data());
  get.add_http_header(Config::query().api_key());
  get.perform();
  HttpStatus http_stat(get);
  http_stat.check(200);

  stts_->set(get.body());
/*
  // The code below attempted to analyse the status data to find out the 
  // progress of the solver. This is an undocumented use and does not seem to 
  // work so far. Achieved here for potential use in the future.
  
  DEB_line(2, stts_);

  const auto& details = stts_->ptree().get_child_optional("details");
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
  get.add_http_header(Config::query().api_key());
  get.perform();
  HttpStatus http_stat(get);
  http_stat.check(200);

  JsonTokens log(get.body());

  bool got_time_data = false;
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

      got_time_data = true;
      const int sol_nmbr = atoi(msg.data() + sol_str_idx + sol_str_len);
      //DEB_line(1, "# solutions found so far: " << sol_nmbr);
      if (sol_nmbr > sol_nmbr_) // new solution(s) found
      {// update the number of solutions and the time of the last solution found
        sol_nmbr_ = sol_nmbr;
        sol_sec_nmbr_ = sec_nmbr;
      }
      stld_sec_nmbr_ = sec_nmbr - sol_sec_nmbr_;
      stop_wtch_.restart();
    }
    log_seq_idx_ = log_item.second.get_child("seqid").get_value<int>() + 1; 
  }
  if (!got_time_data)
  {
    DEB_warning(2, "DOCloud did not provide time and solutions number-"
                   "using internal time counter.");
    stld_sec_nmbr_ += stop_wtch_.restart() / 1000;
  }
}

bool Job::active() const
{
  std::string exct_stts;
  stts_->find_value("executionStatus", exct_stts);
  
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

bool Job::stalled() const 
{ 
  // exit quick if we have a solution, or wait if we don't have one
  const auto& config = Config::query();
  return (sol_nmbr_ > 0 && stld_sec_nmbr_ >= config.feasible_timeout()) 
    || (sol_nmbr_ == 0 && stld_sec_nmbr_ >= config.infeasible_timeout());
}

void Job::abort()
{
  std::string exct_stts;
  stts_->find_value("executionStatus", exct_stts);
  if (exct_stts != "RUNNING")
    return; // already aborted or aborting

  cURLpp::Delete del;
  THROW_OUTCOME_if(!del.valid(), TODO); //Failed to initialize the request
  const std::string url = url_ + "/execute";
  del.set_url(url.data());
  del.add_http_header(Config::query().api_key());
  del.perform();

  HttpStatus http_stat(del);
  http_stat.check(204);
}

void Job::wait()
{
  do 
  {
    std::this_thread::sleep_for(std::chrono::seconds(1));
    sync_status();
    sync_log();
    if (stalled())
      abort(); // The waiting loop must continue until the timer request has been
               // processed by the server.
  } while (active());
}

double Job::solution(std::vector<double>& _x) const
{
  DEB_enter_func;

  // check the solution status (assume it's synchronized already)

  // What are the possible values for this??
  std::string slv_stts;
  stts_->find_value("solveStatus", slv_stts);
  DEB_line(2, "solveStatus=" << slv_stts);
  
  cURLpp::Get get;
  THROW_OUTCOME_if(!get.valid(), TODO); //Failed to initialize the request

  auto url = url_ + "/attachments/solution.json/blob";
  get.set_url(url.data());
  get.add_http_header(Config::query().api_key());
  get.perform();

  HttpStatus http_stat(get);
  if (http_stat.code() == 404 &&
      get.body().find("CIVOC5102E") != std::string::npos)
  {
    // The mixed integer optimization has not found any solution.
    _x.clear();
    return 0;
  }
  http_stat.check(200);

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
  const auto obj_val = 
    bdy_tkns.ptree().get_child("CPLEXSolution.header.objectiveValue").get_value<double>();

  DEB_line(3, "X=" << _x);
  return obj_val;
}

} // namespace DOcloud 
} // namespace COMISO 

#endif // COMISO_DOCLOUD_AVAILABLE
