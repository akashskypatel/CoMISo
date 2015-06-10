#ifndef COMISO_DOCloudCache_HH
#define COMISO_DOCloudCache_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_DOCLOUD_AVAILABLE

//== INCLUDES =================================================================

#include <string>
#include <vector>

//== NAMESPACES ===============================================================
namespace COMISO {
namespace DOcloud {

// given a .lp file, checks if we have saved the result for the same problem.
// If so they are returned, so the calling function can avoid to compute them.
//

class Cache
{
public:
  Cache() : found_(false) {}

  bool restore_result(
    std::string& _file_name,  // .lp file defining the optimization problem
    std::vector<double>& _x,  // result.
    double& _obj_val);        // objective function value.

  // We can store the result for the given .lp file. This makes sense only we have
  // not found cache data for the given .lp file, and in order to avoid data
  // corruption this function fails if the data have been found.
  void store_result(const std::vector<double>& _x, const double& _obj_val);

  const std::string& get_lp_content() { return lp_file_cnts_; }
private:
  std::string key_;          // String generated from .lp file content data. 
                             // This is a sort of hash key generate form the .lp
                             // file content.
  std::string last_filename_;// Last name we have tried to get cached data. 
  std::string lp_file_cnts_; // Content of the input .lp file.
  bool        found_;        // Remembers if we have found a cache for the input
                             // .lp file.
};

} // namespace DOcloud
} // namespace COMISO
//=============================================================================

#endif // COMISO_DOCLOUD_AVAILABLE
//=============================================================================
#endif // COMISO_DOCloudCache_HH
//=============================================================================
