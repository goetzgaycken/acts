#pragma read                           \
  sourceClass="std::vector<int>"      \
  source=""                           \
  version="[1-]"                      \
  targetClass="std::vector<uint64_t>" \
  target=""                           \
  include="vector,cstdint"

#pragma read                           \
  sourceClass="std::vector< std::vector<std::vector<int> > >"       \
  source=""                                                         \
  version="[1-]"                                                    \
  targetClass="std::vector< std::vector< std::vector<uint64_t> > >" \
  target=""                                                         \
  include="vector,cstdint"
