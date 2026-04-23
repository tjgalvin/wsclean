#ifndef WSCLEAN_MSPROVIDERS_WRITE_TO_MS_H_
#define WSCLEAN_MSPROVIDERS_WRITE_TO_MS_H_

namespace schaapcommon::reordering {
class HandleData;
}

namespace wsclean {

void WriteReorderedToMS(const schaapcommon::reordering::HandleData& handle);

};

#endif
