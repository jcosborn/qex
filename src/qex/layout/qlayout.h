#include "qmp.h"

typedef struct llist {
  void *value;
  struct llist *next;
} llist;

// k = sum_i sum_j ((x[i]/d[i][j])%m[i][j])*f[i][j]
// x[i] = sum_j ((k/f[i][j])%m[i][j])*d[i][j]
// parity?

typedef struct {
  int nDim;
  int *physGeom;
  int *rankGeom;
  int *innerGeom; //wrap
  int *outerGeom; //wls
  int *localGeom;
  int physVol;
  int nEven;
  int nOdd;
  int nSites;
  int nEvenOuter;
  int nOddOuter;
  int nSitesOuter;
  int nSitesInner;
  int innerCb;
  int innerCbDir;
  llist *shifts;
  int nranks;
  int myrank;
} Layout;

typedef struct {
  int rank;
  int index;
} LayoutIndex;

typedef struct {
  int begin;
  int end;
  int beginOuter;
  int endOuter;
} Subset;

void layoutSetup(Layout *l);
void layoutIndex(Layout *l, LayoutIndex *li, int coords[]);
void layoutCoord(Layout *l, int coords[], LayoutIndex *li);
void layoutShift(Layout *l, LayoutIndex *li, LayoutIndex *li2, int disp[]);
void layoutSubset(Subset *s, Layout *l, char *sub);

typedef struct {
  int myRank;
  int nIndices;
  int *srcRanks;
  int *srcIndices;
  int nRecvDests;
  int nSendIndices;
  int *sendSrcIndices;
  int *sendDestRanks;
  int *sendDestIndices;
} GatherDescription;

typedef struct {
  GatherDescription *gd;
  int myRank;
  int nIndices;
  int *srcIndices;
  int nRecvRanks;
  int *recvRanks;
  int *recvRankSizes;
  int *recvRankOffsets;
  int recvSize;
  int nRecvDests;
  int *recvDestIndices;
  int *recvBufIndices;
  int nSendRanks;
  int *sendRanks;
  int *sendRankSizes;
  int *sendRankOffsets;
  int sendSize; // same as nSendIndices
  int nSendIndices;
  int *sendIndices;
} GatherIndices;

// per gather:
//  pidx
//  recv
// combined:
//  send*

typedef struct {
  GatherIndices *gi;
  int *disp;
  int *sidx;
  int *pidx;
  int nRecvRanks;
  int *recvRanks;
  int *recvRankSizes;
  int *recvRankSizes1;
  int *recvRankOffsets;
  int *recvRankOffsets1;
  int nRecvSites;
  int nRecvSites1;
  int nRecvDests;
  int *recvDests;
  int *recvLocalSrcs;
  int *recvRemoteSrcs;
  int nSendRanks;
  int *sendRanks;
  int *sendRankSizes;
  int *sendRankSizes1;
  int *sendRankOffsets;
  int *sendRankOffsets1;
  int nSendSites;
  int nSendSites1;
  int *sendSites;
  int vv;
  //int offr, lenr, nthreads;
  int perm;
  int pack;
  int blend;
  //QMP_msgmem_t sqmpmem;
  //QMP_msghandle_t smsg;
  //QMP_msgmem_t rqmpmem;
  //QMP_msghandle_t rmsg;
  //QMP_msghandle_t pairmsg;
} ShiftIndices;

typedef struct {
  QMP_msgmem_t sqmpmem;
  QMP_msghandle_t smsg;
  QMP_msgmem_t rqmpmem;
  QMP_msghandle_t rmsg;
  QMP_msghandle_t pairmsg;
  char *sbuf;
  char *rbuf;
  int sbufSize;
  int rbufSize;
  int first;
  int *offr;
  int *lenr;
  int *nthreads;
} ShiftBuf;

typedef void GatherMap(int *srcRank, int *srcIdx, int dstRank, int *dstIdx, void *args);
