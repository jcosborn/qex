#define _POSIX_C_SOURCE 200112L
#include <stdlib.h>
#include <stdio.h>
#include "qlayout.h"

#define PAIR

#define MAXTHREADS 512
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

static void *
myalloc(size_t size)
{
  size_t align = 64;
#if 0
  char *a = malloc(size+32);
  int o = ((size_t)a) & 31;
  o = (32 - o) & 31;
  char *b = a + o;
  return (void *)b;
#else
  void *b = NULL;
  int err = posix_memalign(&b, align, size);
  if(err) {
    printf("posix_memalign(&%p, %li, %li) failed!\n", b, align, size);
    exit(1);
  }
  return b;
#endif
}

void
prepareShiftBufs(ShiftBuf *sb[], ShiftIndices *si[], int n, int esize)
{
  int sbs=0, rbs=0;
  for(int i=0; i<n; i++) {
    sbs += esize*si[i]->nSendSites1;
    rbs += esize*si[i]->nRecvSites1;
  }
  for(int i=0; i<n; i++) {
    sb[i]->sbufSize = sbs;
    sb[i]->rbufSize = rbs;
    sb[i]->smsg = NULL;
    sb[i]->rmsg = NULL;
    sb[i]->first = 0;
    sb[i]->offr = myalloc(MAXTHREADS*sizeof(int));
    sb[i]->lenr = myalloc(MAXTHREADS*sizeof(int));
    sb[i]->nthreads = myalloc(MAXTHREADS*sizeof(int));
    for(int j=0; j<MAXTHREADS; j++) sb[i]->nthreads[j] = 0;
  }
  sb[0]->first = 1;
  if(sbs>0) {
    void *sbuf = myalloc(sbs);
    //printf("sbuf: %p\n", sbuf);
    for(int i=0; i<n; i++) {
      sb[i]->sbuf = sbuf;
      if(si[i]->nSendRanks>0) {
	sb[i]->sqmpmem = QMP_declare_msgmem(sbuf, sbs);
	sb[i]->smsg =
	  QMP_declare_send_to(sb[i]->sqmpmem, si[i]->sendRanks[0], 0);
	//printf("send: to: %i\tsize: %i\t%p\t%p\n",si[i]->sendRanks[0],sbs,sbuf,sb[i]->smsg);
	//fflush(stdout);
      }
    }
  }
  //printf("rbs: %i\n", rbs);
  if(rbs>0) {
    void *rbuf = myalloc(rbs);
    for(int i=0; i<n; i++) {
      sb[i]->rbuf = rbuf;
      if(si[i]->nRecvRanks>0) {
	sb[i]->rqmpmem = QMP_declare_msgmem(rbuf, rbs);
	sb[i]->rmsg =
	  QMP_declare_receive_from(sb[i]->rqmpmem, si[i]->recvRanks[0], 0);
	//printf("recv: fr: %i\tsize: %i\t%p\t%p\n",si[i]->recvRanks[0],rbs,rbuf,sb[i]->rmsg);
	//fflush(stdout);
      }
    }
  }
  //printf("startSend: %i to: %i\tsize: %i\n", myrank, si->sendRanks[0], size);
  //fflush(stdout);
#ifdef PAIR
  QMP_msghandle_t p[2*n];
  int nn = 0;
  for(int i=0; i<n; i++) {
    if(sb[i]->rmsg) {
      p[nn] = sb[i]->rmsg;
      nn++;
    }
    if(sb[i]->smsg) {
      p[nn] = sb[i]->smsg;
      nn++;
    }
  }
  QMP_msghandle_t pairmsg = NULL;
  if(nn>0) {
    pairmsg = QMP_declare_send_recv_pairs(p, nn);
  }
  for(int i=0; i<n; i++) {
    sb[i]->pairmsg = pairmsg;
    //printf("pair[%i]: %p\t%p\t%p\n",i,sb[i]->rmsg,sb[i]->smsg,sb[i]->pairmsg);
  }
  //fflush(stdout);
#endif
}

void
prepareShiftBuf(ShiftBuf *sb, ShiftIndices *si, int esize)
{
  prepareShiftBufs(&sb, &si, 1, esize);
}

void
prepareShiftBuf2(ShiftBuf *sb, ShiftIndices *si, int esize)
{
  sb->sbufSize = esize*si->nSendSites1;
  if(sb->sbufSize>0) {
    sb->sbuf = myalloc(sb->sbufSize);
    sb->sqmpmem = QMP_declare_msgmem(sb->sbuf, sb->sbufSize);
    sb->smsg = QMP_declare_send_to(sb->sqmpmem, si->sendRanks[0], 0);
  } else {
    sb->smsg = NULL;
  }

  sb->rbufSize = esize*si->nRecvSites1;
  if(sb->rbufSize>0) {
    sb->rbuf = myalloc(sb->rbufSize);
    sb->rqmpmem = QMP_declare_msgmem(sb->rbuf, sb->rbufSize);
    sb->rmsg = QMP_declare_receive_from(sb->rqmpmem, si->recvRanks[0], 0);
  } else {
    sb->rmsg = NULL;
  }

  //printf("startSend: %i to: %i\tsize: %i\n", myrank, si->sendRanks[0], size);
  //fflush(stdout);
#ifdef PAIR
  int n = 0;
  QMP_msghandle_t p[2];
  if(sb->rbufSize>0) {
    p[n] = sb->rmsg;
    n++;
  }
  if(sb->sbufSize>0) {
    p[n] = sb->smsg;
    n++;
  }
  if(n>0) {
    sb->pairmsg = QMP_declare_send_recv_pairs(p, n);
  } else {
    sb->pairmsg = NULL;
  }
#endif
}

void
startSendBuf(ShiftBuf *sb)
{
  //printf("send: %g\n",*(float *)(sb->sbuf));
#ifdef PAIR
  if(sb->pairmsg) QMP_start(sb->pairmsg);
#else
  if(sb->smsg) QMP_start(sb->smsg);
#endif
}
void
startRecvBuf(ShiftBuf *sb)
{
#ifdef PAIR
#else
  if(sb->rmsg) QMP_start(sb->rmsg);
#endif
}
void
waitSendBuf(ShiftBuf *sb)
{
#ifdef PAIR
#else
  if(sb->smsg) QMP_wait(sb->smsg);
#endif
}
void
waitRecvBuf(ShiftBuf *sb)
{
#ifdef PAIR
  if(sb->pairmsg) QMP_wait(sb->pairmsg);
#else
  if(sb->rmsg) QMP_wait(sb->rmsg);
#endif
  //printf("recv: %g\n",*(float *)(sb->rbuf));
}
void
doneRecvBuf(ShiftBuf *sb)
{
#ifdef PAIR
  if(sb->pairmsg) QMP_clear_to_send(sb->pairmsg, QMP_CTS_READY);
#endif
}
void
freeShiftBufs(ShiftBuf *sb[], int n)
{
  for(int i=0; i<n; i++) {
    free(sb[i]->offr);
    free(sb[i]->lenr);
  }
#ifdef PAIR
  for(int i=0; i<n; i++) {
    if(sb[i]->first && sb[i]->pairmsg) {
      QMP_free_msghandle(sb[i]->pairmsg);
    }
    sb[i]->pairmsg = NULL;
    if(sb[i]->smsg) {
      sb[i]->smsg = NULL;
      QMP_free_msgmem(sb[i]->sqmpmem);
      sb[i]->sqmpmem = NULL;
    }
    if(sb[i]->rmsg) {
      sb[i]->rmsg = NULL;
      QMP_free_msgmem(sb[i]->rqmpmem);
      sb[i]->rqmpmem = NULL;
    }
  }
#else
  for(int i=0; i<n; i++) {
    if(sb[i]->smsg) {
      QMP_free_msghandle(sb[i]->smsg);
      sb[i]->smsg = NULL;
      QMP_free_msgmem(sb[i]->sqmpmem);
      sb[i]->sqmpmem = NULL;
    }
    if(sb[i]->rmsg) {
      QMP_free_msghandle(sb[i]->rmsg);
      sb[i]->rmsg = NULL;
      QMP_free_msgmem(sb[i]->rqmpmem);
      sb[i]->rqmpmem = NULL;
    }
  }
#endif
  for(int i=0; i<n; i++) {
    if(sb[i]->first && sb[i]->sbufSize>0) free(sb[i]->sbuf);
    if(sb[i]->first && sb[i]->rbufSize>0) free(sb[i]->rbuf);
  }
}
void
freeShiftBuf(ShiftBuf *sb)
{
  freeShiftBufs(&sb, 1);
}

#if 0
void
startSend(void *buf, int esize, ShiftIndices *si)
{
  if(si->sqmpmem==NULL) {
    int size = esize*si->sendRankSizes1[0];
    //printf("startSend: %i to: %i\tsize: %i\n", myrank, si->sendRanks[0], size);
    //fflush(stdout);
    si->sqmpmem = QMP_declare_msgmem(buf, size);
    si->smsg = QMP_declare_send_to(si->sqmpmem, si->sendRanks[0], 0);
  }
#ifdef PAIR
  if(si->pairmsg==NULL) {
    QMP_msghandle_t p[2] = {si->rmsg,si->smsg};
    si->pairmsg = QMP_declare_send_recv_pairs(p, 2);
  }
  QMP_start(si->pairmsg);
#else
  QMP_start(si->smsg);
#endif
}

void
startRecv(void *buf, int esize, ShiftIndices *si)
{
  if(si->rqmpmem==NULL) {
    int size = esize*si->recvRankSizes1[0];
    //printf("startRecv: %i from: %i\tsize: %i\n", myrank, si->recvRanks[0], size);
    //fflush(stdout);
    si->rqmpmem = QMP_declare_msgmem(buf, size);
    si->rmsg = QMP_declare_receive_from(si->rqmpmem, si->recvRanks[0], 0);
  }
#ifdef PAIR
#else
  QMP_start(si->rmsg);
#endif
}

void
waitSend(ShiftIndices *si)
{
#ifdef PAIR
#else
  QMP_wait(si->smsg);
#endif
}
void
waitRecv(ShiftIndices *si)
{
#ifdef PAIR
  QMP_wait(si->pairmsg);
#else
  QMP_wait(si->rmsg);
#endif
}
void
doneRecvBuf2(ShiftIndices *si)
{
#ifdef PAIR
  QMP_clear_to_send(si->pairmsg, QMP_CTS_READY);
#endif
}

void
freeSend(ShiftIndices *si)
{
#ifdef PAIR
#else
  QMP_free_msghandle(si->smsg);
  QMP_free_msgmem(si->sqmpmem);
  si->sqmpmem = NULL;
#endif
}
void
freeRecv(ShiftIndices *si)
{
#ifdef PAIR
  QMP_free_msghandle(si->pairmsg);
  si->pairmsg = NULL;
  QMP_free_msgmem(si->rqmpmem);
  si->rqmpmem = NULL;
  QMP_free_msgmem(si->sqmpmem);
  si->sqmpmem = NULL;
#else
  QMP_free_msghandle(si->rmsg);
  QMP_free_msgmem(si->rqmpmem);
  si->rqmpmem = NULL;
#endif
}
#endif

typedef struct {
  Layout *l;
  int *disp;
  int parity;
} mapargs;

static void
map(int *sr, int *si, int dr, int *di, void *args)
{
  mapargs *ma = (mapargs *)args;
  Layout *l = ma->l;
  int nd = l->nDim;
  if(*di>=0) {
    int x[nd];
    LayoutIndex dli, sli;
    dli.rank = dr;
    dli.index = *di;
    layoutCoord(l, x, &dli);
    int y[nd];
    int p = 0;
    for(int k=0; k<nd; k++) {
      p += x[k];
      y[k] = (x[k] - ma->disp[k] + l->physGeom[k])%l->physGeom[k];
    }
    if(ma->parity>=0 && (p&1)!=ma->parity) {
      *sr = -1;
      *si = -1;
    } else {
      layoutIndex(l, &sli, y);
      *sr = sli.rank;
      *si = sli.index;
    }
#if 0
    if(myrank==1) {
      printf("%i %i\n", *sr, *si);
      printf("%i %i %i %i -> %i %i %i %i\n",y[0],y[1],y[2],y[3],x[0],x[1],x[2],x[3]);
    }
#endif
  } else {
    // search for site after or including di0 from rank sr to dr
    int di0 = -(*di+1);
    while(di0<l->nSites) {
      int sr0;
      map(&sr0, si, dr, &di0, args);
      if(sr0==*sr) {
	*di = di0;
	return;
      }
      di0++;
    }
    *si = -1;
  }
}

// nRecvRanks (remote ranks)
// start recvs
// nSendRanks
// start sends
// local + perm
// recv buf

// nSendRanks
// - sendRanks
// - nSendPacks
// - - sendPacks
// - - nSendSites
// - - - sendSites

#define SUB2PAR(s) ((s)[0]=='e'?0:((s)[0]=='o'?1:-1))

void
makeGDFromShiftSubs(GatherDescription *gd, Layout *l, int *disps[],
		    char *subs[], int ndisps)
{
  int myRank = l->myrank;
  int myndi = l->nSites;
  int nndi = ndisps*myndi;
  mapargs args;
  args.l = l;

  int *sidx = myalloc(nndi*sizeof(int));
  int *srank = myalloc(nndi*sizeof(int));

  // find shift sources
  int nRecvDests = 0;
  for(int n=0; n<ndisps; n++) {
    int n0 = n*myndi;
    args.disp = disps[n];
    args.parity = SUB2PAR(subs[n]);
#pragma omp parallel for reduction(+:nRecvDests)
    for(int di=0; di<myndi; di++) {
      int sr, si, di0=di;
      map(&sr, &si, myRank, &di0, &args);
      srank[n0+di] = sr;
      sidx[n0+di] = si;
      if(sr != myRank) nRecvDests++;
    }
  }
  gd->myRank = myRank;
  gd->nIndices = nndi;
  gd->srcRanks = srank;
  gd->srcIndices = sidx;
  gd->nRecvDests = nRecvDests;

  // use inverse map
  int nd = l->nDim;
  int dispi[nd];
  args.disp = dispi;
  ARRAY_CREATE(int, sendSrcIndices);
  ARRAY_CREATE(int, sendDestRanks);
  ARRAY_CREATE(int, sendDestIndices);
  int tlen[MAXTHREADS];
  // find who to send to
  for(int n=0; n<ndisps; n++) {
    int n0 = n*myndi;
    int sp = 0;
    for(int i=0; i<nd; i++) {
      sp += abs(disps[n][i]);
      dispi[i] = -disps[n][i];
    }
    args.parity = SUB2PAR(subs[n]);
    if((sp&1)==1 && args.parity>=0) args.parity = 1 - args.parity;
    //#pragma omp parallel
    {
      //int tid = THREADNUM;
      //int nid = NUMTHREADS;
      int tid = 0;
      int nid = 1;
      ARRAY_CREATE(int, sendSrcIndicesT);
      ARRAY_CREATE(int, sendDestRanksT);
      ARRAY_CREATE(int, sendDestIndicesT);
      //#pragma omp for
      for(int di=0; di<myndi; di++) {
	int dr=myRank, sr, si;
	map(&sr, &si, dr, &di, &args);
	if(sr>=0 && si>=0 && sr!=myRank) {
	  if(tid==0) {
	    ARRAY_APPEND(int, sendSrcIndices, di);
	    ARRAY_APPEND(int, sendDestRanks, sr);
	    ARRAY_APPEND(int, sendDestIndices, n0 + si);
	  } else {
	    ARRAY_APPEND(int, sendSrcIndicesT, di);
	    ARRAY_APPEND(int, sendDestRanksT, sr);
	    ARRAY_APPEND(int, sendDestIndicesT, n0 + si);
	  }
	}
      }
      tlen[tid] = nsendSrcIndicesT;
      //TBARRIER;
      int i0 = 0;
      for(int i=0; i<tid; i++) i0 += tlen[i];
      if(tid==nid-1) {
	int len = i0 + nsendSrcIndicesT;
	ARRAY_GROW(int, sendSrcIndices, len);
	ARRAY_GROW(int, sendDestRanks, len);
	ARRAY_GROW(int, sendDestIndices, len);
      }
      i0 += nsendSrcIndices;
      //TBARRIER;
      for(int i=0; i<nsendSrcIndicesT; i++) {
	sendSrcIndices[i0+i] = sendSrcIndicesT[i];
	sendDestRanks[i0+i] = sendDestRanksT[i];
	sendDestIndices[i0+i] = sendDestIndicesT[i];
      }
      if(tid==nid-1) {
	int len = i0 + nsendSrcIndicesT;
	nsendSrcIndices = len;
	nsendDestRanks = len;
	nsendDestIndices = len;
      }
      free(sendSrcIndicesT);
      free(sendDestRanksT);
      free(sendDestIndicesT);
    } // end parallel
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
makeGDFromShifts(GatherDescription *gd, Layout *l, int *disps[], int ndisps)
{
  char *subs[ndisps];
  for(int i=0; i<ndisps; i++) subs[i] = "all";
  makeGDFromShiftSubs(gd, l, disps, subs, ndisps);
}

void
makeShiftMultiSub(ShiftIndices *si[], Layout *l, int *disp[],
		  char *subs[], int ndisp)
{
  int myRank = l->myrank;
  int nd = l->nDim;
  int vvol = l->nSitesOuter;

  GatherIndices *gi = myalloc(sizeof(GatherIndices));
  for(int n=0; n<ndisp; n++) {
    si[n]->gi = gi;
    si[n]->disp = myalloc(nd*sizeof(int));
    for(int i=0; i<nd; i++) {
      si[n]->disp[i] = disp[n][i];
    }
    si[n]->pidx = myalloc(vvol*sizeof(int));
    si[n]->sidx = myalloc(vvol*sizeof(int));
    si[n]->sendSites = myalloc(vvol*sizeof(int));
    for(int i=0; i<vvol; i++) {
      si[n]->pidx[i] = -1;
      si[n]->sidx[i] = -1;
    }
  }

  //mapmargs args;
  //args.l = l;
  //args.disp = disp;
  //args.ndisp = ndisp;
  //makeGather(gi, mapm, &args,l->nranks,l->nranks,l->nSites*ndisp,l->myrank);
  GatherDescription *gd = myalloc(sizeof(GatherDescription));
  makeGDFromShiftSubs(gd, l, disp, subs, ndisp);
  makeGatherFromGD(gi, gd);
#if 0
  if(myrank==0) {
    printf("sidx: %p\n", gi->srcIndices);
    for(int i=0; i<9; i++) {
      printf("%i:\tsidx: %i\n", i, gi->srcIndices[i]);
    }
    fflush(stdout);
    QMP_barrier();
  }
#endif

  //int si0 = 0;
  int si0 = ndisp - 1;
  int vvs=0, perm=0, pack=0;
  //TRACE_ALL;
  if(gi->nSendIndices>0) {
    //if(myrank==1){printf("nss: %i\n", gi->nSendIndices);fflush(stdout);}
    pack = gi->sendIndices[0] % l->nSitesInner;
    if(pack==0) {
      int i=1;
      while( (i<gi->nSendIndices) &&
	     (i<l->nSitesInner) &&
	     (gi->sendIndices[i]==gi->sendIndices[0]+i)) i++;
      pack = -(i % l->nSitesInner);
    }
    int ssi0 = -1;
    for(int i=0; i<gi->nSendIndices; i++) {
      int ss = gi->sendIndices[i];
      int ssi = ss / l->nSitesInner;
      if(ssi!=ssi0) {
	si[si0]->sendSites[vvs] = ssi;
	vvs++;
	ssi0 = ssi;
	if(vvs>vvol) {
	  printf("vvs(%i)>vvol(%i)\n", vvs, vvol);
	  if(myRank==0) {
	    for(int i=0; i<gi->nSendIndices; i++) {
	      printf("%i\t%i\n", i, gi->sendIndices[i]);
	    }
	  }
	  fflush(stdout);
	  QMP_barrier();
	  exit(1);
	}
	//if(myrank==1){printf("vvs: %i\tss: %i\tssi: %i\n",vvs,ss,ssi);fflush(stdout);}
      }
    }
  }
  //return;

#if 0
  if(myrank==1) {
    printf("sidx: %p\n", gi->sidx);fflush(stdout);
    printf("sidx[0]: %i\n", gi->sidx[0]);
  }
#endif
  for(int i=0; i<ndisp; i++) {
    si[i]->nSendRanks = 0;
    si[i]->nSendSites1 = 0;
  }
  si[si0]->nSendRanks = gi->nSendRanks;
  si[si0]->nSendSites = vvs;
  si[si0]->nSendSites1 = gi->nSendIndices;
  if(gi->nSendRanks>0) {
    si[si0]->sendRanks = gi->sendRanks;
    si[si0]->sendRankSizes = myalloc(si[si0]->nSendRanks*sizeof(int));
    si[si0]->sendRankSizes1 = gi->sendRankSizes;
    si[si0]->sendRankOffsets = myalloc(si[si0]->nSendRanks*sizeof(int));
    si[si0]->sendRankOffsets1 = gi->sendRankOffsets;
    si[si0]->sendRankSizes[0] = vvs;
    si[si0]->sendRankOffsets[0] = 0;
  }
  //TRACE_ALL;

#if 0
  if(myrank==1) {
    printf("disp:");
    for(int i=0; i<nd; i++) printf(" %i", disp[i]);
    printf("\n");
    printf("sidx: %p\n", gi->sidx);fflush(stdout);
    printf("sidx[0]: %i\n", gi->sidx[0]);
  }
#endif
  int nrsites=0, nrdests[ndisp];
  for(int i=0; i<ndisp; i++) nrdests[i] = 0;
  for(int i=0; i<vvol*ndisp; i++) {
    //if(myrank==1){printf("%i\n", i);fflush(stdout);}
    int dd = i/l->nSitesOuter;
    int ix = i%l->nSitesOuter;
    int k0 = i*l->nSitesInner;
    int recv = 0;
    int rbi = 0;
    for(int ii=0; ii<l->nSitesInner; ii++) {
      int k = k0 + ii;
      int s = gi->srcIndices[k];
      if(s==-1) { recv = -1; break; }
      if(s<0) {
	recv++;
	if(rbi==0) rbi = s;
      }
    }
    if(recv<0) {
      si[dd]->pidx[ix] = -1;
      si[dd]->sidx[ix] = -1;
    } else if(recv==0) {
      si[dd]->pidx[ix] = gi->srcIndices[k0]/l->nSitesInner;
      si[dd]->sidx[ix] = gi->srcIndices[k0]/l->nSitesInner;
      int p = gi->srcIndices[k0] % l->nSitesInner;
      if(p!=0) {
	perm = p;
	//si->sidx[i] = -vvs-1;
	si[dd]->pidx[ix] = -(si[dd]->pidx[ix])-2;
	//vvs++;
      }
    } else {
      rbi = -(rbi+2);
      rbi = (2*rbi)/l->nSitesInner;
      if(pack==0) rbi /= 2;
      si[dd]->sidx[ix] = -rbi-2;
      //nrsites++;
      nrdests[dd]++;
    }
  }
  //TRACE_ALL;
  nrsites = gi->recvSize/l->nSitesInner;
  if(pack!=0) nrsites *= 2;

  for(int i=0; i<ndisp; i++) {
    si[i]->nRecvRanks = 0;
    si[i]->nRecvSites1 = 0;
  }
  si[0]->nRecvRanks = gi->nRecvRanks;
  si[0]->nRecvSites = nrsites;
  si[0]->nRecvSites1 = gi->recvSize;
  if(gi->nRecvRanks>0) {
    si[0]->recvRanks = gi->recvRanks;
    si[0]->recvRankSizes = myalloc(si[0]->nRecvRanks*sizeof(int));
    si[0]->recvRankSizes1 = gi->recvRankSizes;
    si[0]->recvRankOffsets = myalloc(si[0]->nRecvRanks*sizeof(int));
    si[0]->recvRankOffsets1 = gi->recvRankOffsets;
    si[0]->recvRankSizes[0] = nrsites;
    si[0]->recvRankOffsets[0] = 0;
  }
  //printf("nSendRanks: %i\tnRecvRanks: %i\n", si->nSendRanks, si->nRecvRanks);

  for(int n=0; n<ndisp; n++) {
    si[n]->nRecvDests = nrdests[n];
    if(nrdests[n]>0) {
      si[n]->recvDests = myalloc(nrdests[n]*sizeof(int));
      si[n]->recvLocalSrcs = myalloc(nrdests[n]*sizeof(int));
      si[n]->recvRemoteSrcs = myalloc(nrdests[n]*sizeof(int));
      int j = 0;
      for(int i=0; i<vvol; i++) {
	if(si[n]->sidx[i]<-1) {
	  int k = -(si[n]->sidx[i]+2);
	  si[n]->recvDests[j] = i;
	  si[n]->recvRemoteSrcs[j] = k;
	  si[n]->recvLocalSrcs[j] = 0;
	  for(int i0=0; i0<l->nSitesInner; i0++) {
	    int ii = n*l->nSites + i*l->nSitesInner + i0;
	    int gs = gi->srcIndices[ii];
	    if(gs>=0) {
	      si[n]->recvLocalSrcs[j] = gs/l->nSitesInner;
	      break;
	    }
	  }
	  j++;
	  if(j>nrdests[n]) {printf("j(%i)>nrdests[%i](%i)\n",j,n,nrdests[n]);fflush(stdout);}
	}
      }
    }
    si[n]->vv = vvol;
    si[n]->perm = perm;
    si[n]->pack = pack;
    si[n]->blend = pack;
    //si[n]->offr = 0;
    //si[n]->lenr = 0;
    //si[n]->nthreads = 0;

    //si[n]->sqmpmem = NULL;
    //si[n]->rqmpmem = NULL;
    //si[n]->pairmsg = NULL;

    //printf("%i nsend: %i  nrecv: %i\n", myrank, si[n]->nSendSites, si[n]->nRecvSites);
  }
  //printf("disp:");
  //for(int i=0; i<nd; i++) printf(" %i", disp[i]);
  //printf("\n");
  //printf("  perm: %i\n", perm);
}

void
makeShiftMulti(ShiftIndices *si[], Layout *l, int *disp[], int ndisp)
{
  char *subs[ndisp];
  for(int i=0; i<ndisp; i++) subs[i] = "all";
  makeShiftMultiSub(si, l, disp, subs, ndisp);
}

void
makeShift(ShiftIndices *si, Layout *l, int *disp)
{
  makeShiftMulti(&si, l, &disp, 1);
}

void
makeShiftSub(ShiftIndices *si, Layout *l, int *disp, char *sub)
{
  makeShiftMultiSub(&si, l, &disp, &sub, 1);
}
