# SimTOST 1.0.2

## Bug Fixes
- **Fixed runtime error:** Resolved an issue where a negative value (`-1`) was incorrectly assigned to an `unsigned int`, leading to a runtime error:  
  _"-1 is outside the range of representable values of type 'unsigned int'"_.  
  The fix involved replacing `arma::uvec` with `arma::ivec` to correctly handle signed integers in relevant functions.

# SimTOST 1.0.1
- Fixed CRAN review issues: expanded description, added references, documented function outputs.

# SimTOST 1.0.0
* Initial CRAN submission.

# SimTOST 0.6.0

# SimTOST 0.5.0

# SimTOST 0.4.0

# SimTOST 0.3.0

# SimTOST 0.2.0

# SimTOST 0.1.0


