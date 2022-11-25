#pragma link C++ class std::vector<std::vector<std::vector<int> > >+;
#pragma link C++ class std::vector<std::vector<std::vector<int64_t> > >+;
#pragma link C++ class std::vector<std::vector<std::vector<uint64_t> > >+;

// to use the same storage for athena barcode branches int and Acts barcodes uint64_t
#pragma read                                                      \
  sourceClass="std::vector<std::vector<std::vector<int> > >"      \
  source=""                                                       \
  version="[1-]"                                                  \
  targetClass="std::vector<std::vector<std::vector<uint64_t> > >" \
  target=""                                                       \
  include="vector,cstdint"

#pragma read                          \
  sourceClass="std::vector<int>"      \
  source=""                           \
  version="[1-]"                      \
  targetClass="std::vector<uint64_t>" \
  target=""                           \
  include="vector,cstdint"

