#include <stdlib.h>
#include <stdio.h>
#include "qlayout.h"

#define BEGIN_FATAL (void)0
#define PRINT_FATAL(...) printf0(__VA_ARGS__)
#define END_FATAL exit(-1)
#define CAT(a,b) CATX(a,b)
#define CATX(a,b) a ## b
#define ARRAY_CREATE(tp,nm) int CAT(n,nm)=0, CAT(nm,len)=0; tp *nm=NULL
#define ARRAY_GROW(tp,nm,ln)                    \
  if(CAT(n,nm)+(ln)>CAT(nm,len)) {              \
    CAT(nm,len) += (ln);                        \
    nm = realloc(nm, CAT(nm,len)*sizeof(tp));   \
  }
#define ARRAY_APPEND(tp,nm,vl)                  \
  if(CAT(n,nm)==CAT(nm,len)) {                  \
    CAT(nm,len) *= 2;                           \
    if(CAT(nm,len)==0) CAT(nm,len) = 16;        \
    nm = realloc(nm, CAT(nm,len)*sizeof(tp));   \
  }                                             \
  nm[CAT(n,nm)] = vl;                           \
  CAT(n,nm)++
#define ARRAY_COPY(tp,ds,nm) \
  do { for(int _i=0; _i<CAT(n,nm); _i++) (ds)[_i] = (nm)[_i]; } while(0)
#define ARRAY_CLONE(tp,ds,nm) \
  do { (ds) = myalloc(CAT(n,nm)*sizeof(tp)); ARRAY_COPY(tp,ds,nm); } while(0)
#define printf0(...) if(myRank==0) { printf(__VA_ARGS__); fflush(stdout); }
#define myalloc malloc

// map(&sr,&si,dr>=0,&di>=0) -> sr,si
// map(&sr>=0,&si,dr>=0,&di0<0) ->
//   si>=0 then (sr,si)->(dr,di) and di>=-(di0+1) is smallest such di
//   si<0 if sr doesn't send to dr

// each dest site has zero or one source sites
// each source site can have any number of dest sites
// map(&sr,&si,&dr>=0,&di>=0,dn<0) -> sr,si
// map(&sr>=0,&si>=0,&dr,&di,dn>=0) -> dr,di for dest number dn

// GatherDescription
// pass in myRank, nIndices, srcRanks, srcIndices, 
//  nSendIndices, sendSrcIndices, sendDestRanks, sendDestIndices (not sorted)

// merge (add?)  + (need destIndexOffset per rank/gather)
// compose       *
// filter (on dest sites)

// -> getSendList (1 for multi)
// -> getRecvList (1 for multi)
// -> getDestList (N for multi)
// -> getSrcList (N for multi) (needs recv sites)

// getRecvInfo
// getSendInfo
//  multi

void
mergeGatherDescriptions(GatherDescription *gd, GatherDescription *gds, int n)
{
  int ni=0, nsi=0;
  for(int i=0; i<n; i++) {
    ni += gds[i].nIndices;
    nsi += gds[i].nSendIndices;
    if(gds[i].myRank!=gds[0].myRank) {
      int myRank = 0;
      BEGIN_FATAL;
      PRINT_FATAL("ranks don't match: gds[%i].myRank(%i)!=gds[0].myRank(%i)\n",
		  i, gds[i].myRank, gds[0].myRank);
      END_FATAL;
    }
  }
  int *sr = myalloc(ni*sizeof(int));
  int *si = myalloc(ni*sizeof(int));
  ni = 0;
  for(int i=0; i<n; i++) {
    for(int j=0; j<gds[i].nIndices; j++) {
      sr[ni] = gds[i].srcRanks[j];
      si[ni] = gds[i].srcIndices[j];
      ni++;
    }
  }
  int *ssi = myalloc(nsi*sizeof(int));
  int *sdr = myalloc(nsi*sizeof(int));
  int *sdi = myalloc(nsi*sizeof(int));
  nsi = 0;
  for(int i=0; i<n; i++) {
    for(int j=0; j<gds[i].nSendIndices; j++) {
      ssi[nsi] = gds[i].sendSrcIndices[j];
      sdr[nsi] = gds[i].sendDestRanks[j];
      sdi[nsi] = gds[i].sendDestIndices[j];
      nsi++;
    }
  }
  gd->myRank = gds[0].myRank;
  gd->nIndices = ni;
  gd->srcRanks = sr;
  gd->srcIndices = si;
  gd->nSendIndices = nsi;
  gd->sendSrcIndices = ssi;
  gd->sendDestRanks = sdr;
  gd->sendDestIndices = sdi;
}

static int
cyclicComp(int a, int b, int zero)
{
  int c = abs(a) + abs(b) + 1;
  if(a<zero) a += c;
  if(b<zero) b += c;
  return a-b;
}

static GatherDescription *g_gd;
static int *g_sd;
#if 1
static int
sortSd(const void *a, const void *b)
{
  const int *pa = (const int *)a;
  const int *pb = (const int *)b;
  int mr = g_gd->myRank;
  int ra = pa[0];
  int rb = pb[0];
  int rr = cyclicComp(ra,rb,mr);
  if(rr==0) {
    int ia = pa[2];
    int ib = pb[2];
    rr = ia - ib; 
  }
  return rr;
}
static int
sortSd2(const void *a, const void *b)
{
  const int pa = 3*(*(const int *)a);
  const int pb = 3*(*(const int *)b);
  int mr = g_gd->myRank;
  int ra = g_sd[pa];
  int rb = g_sd[pb];
  int rr = cyclicComp(ra,rb,mr);
  if(rr==0) {
    int ia = g_sd[pa+1];
    int ib = g_sd[pb+1];
    rr = ia - ib; 
    if(rr==0) {
      int ja = g_sd[pa+2];
      int jb = g_sd[pb+2];
      rr = ja - jb; 
    }
  }
  return rr;
}
#endif

// uses:  myRank, nIndices, srcRanks, srcIndices, 
void
makeRecvInfo(GatherIndices *gi, GatherDescription *gd)
{
  int mr = gd->myRank;
  int ni = gd->nIndices;
  gi->gd = gd;
  gi->myRank = mr;
  gi->nIndices = ni;
  gi->srcIndices = NULL;
  gi->nRecvRanks = 0;
  gi->recvRanks = NULL;
  gi->recvRankSizes = NULL;
  gi->recvRankOffsets = NULL;
  gi->recvSize = 0;
  gi->nRecvDests = 0;
  gi->recvDestIndices = NULL;
  gi->recvBufIndices = NULL;
  if(ni==0) return;

#if 1
  // sort: srcRank, destIndex
  // perm(p): srcRank, srcIndex, destIndex  : a[p[i]] < a[p[j]] (i<j)
  // add firstDestIndex[p[i]] = p[i-?]
  // rdi[i] = destIndex[i]
  // if(firstDestIndex
  // rbi[i] = j
  // sd: srcRanks, srcIndices, destIndices
  int *sd = myalloc(3*gd->nRecvDests*sizeof(int));
  int *si = myalloc(ni*sizeof(int));
  int nsd = 0;
  //#pragma omp parallel for
  for(int i=0; i<ni; i++) {
    int r = gd->srcRanks[i];
    if(r<0) {
      si[i] = -1;
    } else if(r==mr) { //local
      si[i] = gd->srcIndices[i];
    } else { // remote
      //#pragma omp critical
      {
	sd[nsd] = r;
	sd[nsd+1] = gd->srcIndices[i];
	sd[nsd+2] = i;
	nsd += 3;
      }
    }
  }
  nsd /= 3;
  g_gd = gd;
  qsort(sd, nsd, 3*sizeof(int), sortSd);

  int *p = myalloc(nsd*sizeof(int));
  for(int i=0; i<nsd; i++) p[i] = i;
  g_sd = sd;
  qsort(p, nsd, sizeof(int), sortSd2);

  ARRAY_CREATE(int, recvRanks);
  ARRAY_CREATE(int, recvRankCounts);
  // check for duplicate source indices
  {
    int rr = -1;
    int rrc = 0;
    int k0 = -1;
    int r0 = -1;
    int i0 = -1;
    for(int i=0; i<nsd; i++) {
      int k = p[i];
      int k3 = 3*k;
      int sri = sd[k3];
      int sii = sd[k3+1];
      if(sri!=r0 || sii!=i0) {
	k0 = k;
	r0 = sri;
	i0 = sii;
	if(r0!=rr) {
	  if(i>0) { ARRAY_APPEND(int, recvRankCounts, rrc); }
	  ARRAY_APPEND(int, recvRanks, r0);
	  rr = r0;
	  rrc = 0;
	}
	rrc++;
      }
      sd[k3+1] = k0;
    }
    ARRAY_APPEND(int, recvRankCounts, rrc);
  }

  int rsize = 0;
  int *rr = myalloc(nrecvRanks*sizeof(int));
  int *rrs = myalloc(nrecvRanks*sizeof(int));
  int *rro = myalloc(nrecvRanks*sizeof(int));
  for(int i=0; i<nrecvRanks; i++) {
    rr[i] = recvRanks[i];
    rrs[i] = recvRankCounts[i];
    rro[i] = rsize;
    rsize += recvRankCounts[i];
  }

  // finish si, rbuf
  int *rdi = myalloc(nsd*sizeof(int));
  int *rbi = myalloc(nsd*sizeof(int));
  //#pragma omp parallel for
  int j = 0;
  for(int i=0; i<nsd; i++) {
    int i3 = 3*i;
    int di = sd[i3+2];
    int fdi = sd[i3+1];
    int bi = j;
    if(fdi==i) j++;
    else bi = sd[3*fdi+1];
    sd[i3+1] = bi;
    rdi[i] = di;
    rbi[i] = bi;
    si[di] = -bi-2;
  }

  gi->srcIndices = si;
  gi->nRecvRanks = nrecvRanks;
  gi->recvRanks = rr;
  gi->recvRankSizes = rrs;
  gi->recvRankOffsets = rro;
  gi->recvSize = rsize;
  gi->nRecvDests = nsd;
  gi->recvDestIndices = rdi;
  gi->recvBufIndices = rbi;

  free(sd);
  free(p);
  free(recvRanks);
  free(recvRankCounts);

#else
  ARRAY_CREATE(int, recvRanks);
  ARRAY_CREATE(int, recvRankCounts);
  ARRAY_CREATE(int, rbufSrcRankIndices);
  ARRAY_CREATE(int, rbufSrcIndices);
  ARRAY_CREATE(int, rbufRankIndices);
  ARRAY_CREATE(int, recvDestIndices);
  ARRAY_CREATE(int, recvBufIndices);
  int *si = myalloc(ni*sizeof(int));
  for(int i=0; i<ni; i++) {
    int r = gd->srcRanks[i];
    if(r<0) {
      si[i] = -1;
    } else if(r==mr) { //local
      si[i] = gd->srcIndices[i];
    } else { // remote
      // check if new rank
      int ri = 0;
      while(ri<nrecvRanks && r!=recvRanks[ri]) ri++;
      if(ri==nrecvRanks) {
	ARRAY_APPEND(int, recvRanks, r);
	ARRAY_APPEND(int, recvRankCounts, 0);
      }
      // check if in rbuf
      int sii = gd->srcIndices[i];
      int j = 0;
      while(j<nrbufSrcRankIndices &&
	    (rbufSrcRankIndices[j]!=ri || rbufSrcIndices[j]!=sii)) j++;
      // if not found, add to rbuf
      if(j==nrbufSrcRankIndices) {
	ARRAY_APPEND(int, rbufSrcRankIndices, ri);
	ARRAY_APPEND(int, rbufSrcIndices, sii);
	ARRAY_APPEND(int, rbufRankIndices, recvRankCounts[ri]);
	recvRankCounts[ri]++;
      }
      ARRAY_APPEND(int, recvDestIndices, i);
      //ARRAY_APPEND(int, recvBufIndices, rbufRankIndices[j]);
      ARRAY_APPEND(int, recvBufIndices, j);
    }
  }

  // sort ranks with indexing array
  int *p = myalloc(nrecvRanks*sizeof(int));
  for(int i=0; i<nrecvRanks; i++) p[i] = i;
  for(int i=0; i<nrecvRanks; i++) {
    int ri = i;
    int rv = recvRanks[p[i]];
    for(int j=i+1; j<nrecvRanks; j++) {
      if(cyclicComp(recvRanks[p[j]],rv,mr)<=0) {
	ri = j;
	rv = recvRanks[p[j]];
      }
    }
    int l = p[i];
    p[i] = p[ri];
    p[ri] = l;
  }
  int *pinv = myalloc(nrecvRanks*sizeof(int));
  for(int i=0; i<nrecvRanks; i++) {
    pinv[p[i]] = i;
  }

  int rsize = 0;
  int *rr = myalloc(nrecvRanks*sizeof(int));
  int *rrs = myalloc(nrecvRanks*sizeof(int));
  int *rro = myalloc(nrecvRanks*sizeof(int));
  for(int i=0; i<nrecvRanks; i++) {
    int ri = p[i];
    rr[i] = recvRanks[ri];
    rrs[i] = recvRankCounts[ri];
    rro[i] = rsize;
    rsize += recvRankCounts[ri];
  }

  // finish si, rbuf
  int *rdi = myalloc(nrecvDestIndices*sizeof(int));
  int *rbi = myalloc(nrecvDestIndices*sizeof(int));
#pragma omp parallel for
  for(int i=0; i<nrecvDestIndices; i++) {
    int di = recvDestIndices[i];
    int j = recvBufIndices[i];
    int ri = rbufSrcRankIndices[j];
    int bi = rbufRankIndices[j];
    rdi[i] = di;
    rbi[i] = rro[pinv[ri]] + bi;
    si[di] = -rbi[i]-2;
  }

  //gi->gd = gd;
  //gi->myRank = mr;
  //gi->nIndices = ni;
  gi->srcIndices = si;
  gi->nRecvRanks = nrecvRanks;
  gi->recvRanks = rr;
  gi->recvRankSizes = rrs;
  gi->recvRankOffsets = rro;
  gi->recvSize = rsize;
  gi->nRecvDests = nrecvDestIndices;
  gi->recvDestIndices = rdi;
  gi->recvBufIndices = rbi;

  free(recvRanks);
  free(recvRankCounts);
  free(rbufSrcRankIndices);
  free(rbufSrcIndices);
  free(rbufRankIndices);
  free(recvDestIndices);
  free(recvBufIndices);
  free(p);
  free(pinv);
#endif
}

static int
sortsendSrc(const void *a, const void *b)
{
  const int pa = *(const int *)a;
  const int pb = *(const int *)b;
  int mr = g_gd->myRank;
  int ra = g_gd->sendDestRanks[pa];
  int rb = g_gd->sendDestRanks[pb];
  int rr = cyclicComp(ra,rb,mr);
  if(rr==0) {
    int ia = g_gd->sendSrcIndices[pa];
    int ib = g_gd->sendSrcIndices[pb];
    rr = ia - ib; 
    if(rr==0) {
      int ja = g_gd->sendDestIndices[pa];
      int jb = g_gd->sendDestIndices[pb];
      rr = ja - jb; 
    }
  }
  return rr;
}

static int
sortsend(const void *a, const void *b)
{
  const int pa = *(const int *)a;
  const int pb = *(const int *)b;
  int mr = g_gd->myRank;
  int ra = g_gd->sendDestRanks[pa];
  int rb = g_gd->sendDestRanks[pb];
  int rr = cyclicComp(ra,rb,mr);
  if(rr==0) {
    int ia = g_gd->sendDestIndices[pa];
    int ib = g_gd->sendDestIndices[pb];
    rr = ia - ib; 
  }
  return rr;
}

// uses: myRank, nSendIndices, sendSrcIndices, sendDestRanks, sendDestIndices
void
makeSendInfo(GatherIndices *gi, GatherDescription *gd)
{
  gi->nSendRanks = 0;
  gi->sendRanks = NULL;
  gi->sendRankSizes = NULL;
  gi->sendRankOffsets = NULL;
  gi->sendSize = 0;
  gi->nSendIndices = 0;
  gi->sendIndices = NULL;
  int n = gd->nSendIndices;
  if(n==0) return;

#if 1
  int *p = myalloc(n*sizeof(int));
  for(int i=0; i<n; i++) p[i] = i;
  g_gd = gd;
  qsort(p, n, sizeof(int), sortsendSrc);
  ARRAY_CREATE(int, sendRanks);
  ARRAY_CREATE(int, sendRankCounts);
  int ndup = 0;
  {
    int sr = -1;
    int src = 0;
    int r0 = -1;
    int i0 = -1;
    for(int i=0; i<n; i++) {
      int pi = p[i];
      int sdr = gd->sendDestRanks[pi];
      int ssi = gd->sendSrcIndices[pi];
      if(sdr!=r0 || ssi!=i0) {
	r0 = sdr;
	i0 = ssi;
	if(r0!=sr) {
	  if(i>0) { ARRAY_APPEND(int, sendRankCounts, src); }
	  ARRAY_APPEND(int, sendRanks, r0);
	  sr = r0;
	  src = 0;
	}
	src++;
      } else {
	p[i] = p[ndup];
	p[ndup] = -1;
	ndup++;
      }
    }
    ARRAY_APPEND(int, sendRankCounts, src);
  }

  int ssize = 0;
  int *sr = myalloc(nsendRanks*sizeof(int));
  int *srs = myalloc(nsendRanks*sizeof(int));
  int *sro = myalloc(nsendRanks*sizeof(int));
  for(int i=0; i<nsendRanks; i++) {
    sr[i] = sendRanks[i];
    srs[i] = sendRankCounts[i];
    sro[i] = ssize;
    ssize += sendRankCounts[i];
  }

  int nsend = n - ndup;
  qsort(p+ndup, nsend, sizeof(int), sortsend);

  int *si = myalloc(nsend*sizeof(int));
  for(int i=0; i<nsend; i++) {
    int k = p[ndup+i];
    si[i] = gd->sendSrcIndices[k];
  }

  gi->nSendRanks = nsendRanks;
  gi->sendRanks = sr;
  gi->sendRankSizes = srs;
  gi->sendRankOffsets = sro;
  gi->sendSize = nsend;
  gi->nSendIndices = nsend;
  gi->sendIndices = si;

  free(p);
  free(sendRanks);
  free(sendRankCounts);

#else
  int *p = myalloc(n*sizeof(int));
  for(int i=0; i<n; i++) p[i] = i;
  int mr = gd->myRank;
  for(int i=0; i<n; i++) {
    int k = i;
    int sdr = gd->sendDestRanks[p[i]];
    int sdi = gd->sendDestIndices[p[i]];
    for(int j=i+1; j<n; j++) {
      int rj = gd->sendDestRanks[p[j]];
      int ij = gd->sendDestIndices[p[j]];
      if(((rj==sdr)&&(ij<sdi)) || (cyclicComp(rj,sdr,mr)<0)) {
	k = j;
	sdr = rj;
	sdi = ij;
      }
    }
    int l = p[i];
    p[i] = p[k];
    p[k] = l;
  }

  ARRAY_CREATE(int, sendRanks);
  ARRAY_CREATE(int, sendRankSizes);
  ARRAY_CREATE(int, sendRankOffsets);
  ARRAY_CREATE(int, sendIndices);
  ARRAY_APPEND(int, sendRanks, gd->sendDestRanks[p[0]]);
  ARRAY_APPEND(int, sendRankOffsets, 0);
  for(int i=0; i<n; i++) {
    int r = gd->sendDestRanks[p[i]];
    if(r!=sendRanks[nsendRanks-1]) {
      ARRAY_APPEND(int, sendRankSizes, nsendIndices-sendRankOffsets[nsendRankOffsets-1]);
      ARRAY_APPEND(int, sendRanks, r);
      ARRAY_APPEND(int, sendRankOffsets, nsendIndices);
    }
    int ssi = gd->sendSrcIndices[p[i]];
    int j = sendRankOffsets[nsendRankOffsets-1];
    while(j<nsendIndices && ssi!=sendIndices[j]) j++;
    if(j==nsendIndices) {
      ARRAY_APPEND(int, sendIndices, ssi);
    }
  }
  ARRAY_APPEND(int, sendRankSizes, nsendIndices-sendRankOffsets[nsendRankOffsets-1]);

  gi->nSendRanks = nsendRanks;
  gi->sendRanks = myalloc(nsendRanks*sizeof(int));
  ARRAY_CLONE(int, gi->sendRanks, sendRanks);
  ARRAY_CLONE(int, gi->sendRankSizes, sendRankSizes);
  ARRAY_CLONE(int, gi->sendRankOffsets, sendRankOffsets);
  gi->sendSize = nsendIndices;
  gi->nSendIndices = nsendIndices;
  ARRAY_CLONE(int, gi->sendIndices, sendIndices);

  free(p);
  free(sendRanks);
  free(sendRankSizes);
  free(sendRankOffsets);
  free(sendIndices);
#endif
}

void
makeGD(GatherDescription *gd, GatherMap *map, void *args,
       int nSrcRanks, int nDstRanks, int myndi, int myRank)
{
  int *sidx = myalloc(myndi*sizeof(int));
  int *srank = myalloc(myndi*sizeof(int));
  // find shift sources
  for(int di=0; di<myndi; di++) {
    int sr, si, di0=di;
    map(&sr, &si, myRank, &di0, args);
    srank[di] = sr;
    sidx[di] = si;
  }

  gd->myRank = myRank;
  gd->nIndices = myndi;
  gd->srcRanks = srank;
  gd->srcIndices = sidx;

  ARRAY_CREATE(int, sendSrcIndices);
  ARRAY_CREATE(int, sendDestRanks);
  ARRAY_CREATE(int, sendDestIndices);
  // find who to send to
  for(int dr=0; dr<nDstRanks; dr++) {
    if(dr==myRank) continue;
    int sr=myRank, si, di=-1;
    map(&sr, &si, dr, &di, args);
    while(si>=0) {
      ARRAY_APPEND(int, sendSrcIndices, si);
      ARRAY_APPEND(int, sendDestRanks, dr);
      ARRAY_APPEND(int, sendDestIndices, di);
      di = -di-2;
      map(&sr, &si, dr, &di, args);
    }
  }

  gd->nSendIndices = nsendSrcIndices;
  ARRAY_CLONE(int, gd->sendSrcIndices, sendSrcIndices);
  ARRAY_CLONE(int, gd->sendDestRanks, sendDestRanks);
  ARRAY_CLONE(int, gd->sendDestIndices, sendDestIndices);
  free(sendSrcIndices);
  free(sendDestRanks);
  free(sendDestIndices);
}

void
makeGathersFromGDs(GatherIndices *gi[], GatherDescription *gd[], int n)
{
  for(int i=0; i<n; i++) {
    makeRecvInfo(gi[i], gd[i]);
    makeSendInfo(gi[i], gd[i]);
  }
}

void
makeGatherFromGD(GatherIndices *gi, GatherDescription *gd)
{
  makeRecvInfo(gi, gd);
  makeSendInfo(gi, gd);
}

void
makeGather(GatherIndices *gi, GatherMap *map, void *args,
	   int nSrcRanks, int nDstRanks, int myndi, int myRank)
{
  GatherDescription *gd = myalloc(sizeof(GatherDescription));
  makeGD(gd, map, args, nSrcRanks, nDstRanks, myndi, myRank);
  makeRecvInfo(gi, gd);
  makeSendInfo(gi, gd);
}
