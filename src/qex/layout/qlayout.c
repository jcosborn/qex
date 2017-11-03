#include <stdlib.h>
#include <stdio.h>
#include "qlayout.h"

#define myalloc malloc
#define PRINTV(s,f,v,n) do { printf(s);                 \
    for(int _i=0; _i<n; _i++) printf(" "f, (v)[_i]);    \
    printf("\n"); } while(0)

void
layoutSubset(Subset *s, Layout *l, char *sub)
{
  s->begin = 0;
  s->end = l->nSites;
  s->beginOuter = 0;
  s->endOuter = l->nSitesOuter;
  if(sub[0]=='e') {
    s->end = l->nEven;
    s->endOuter = l->nEvenOuter;
  } else if(sub[0]=='o') {
    s->begin = l->nOdd;
    s->beginOuter = l->nOddOuter;
  }
}

void
layoutSetup(Layout *l)
{
  int nd = l->nDim;
  l->outerGeom = myalloc(nd*sizeof(int));
  l->localGeom = myalloc(nd*sizeof(int));
  int pvol=1, lvol=1, ovol=1, icb=0, icbd=-1;
  for(int i=0; i<nd; i++) {
    l->localGeom[i] = l->physGeom[i]/l->rankGeom[i];
    l->outerGeom[i] = l->localGeom[i]/l->innerGeom[i];
    pvol *= l->physGeom[i];
    lvol *= l->localGeom[i];
    ovol *= l->outerGeom[i];
    if(l->innerGeom[i]>1 && (l->outerGeom[i]&1)==1) icb++;
    if(l->innerGeom[i]==1 && (l->outerGeom[i]&1)==0) icbd = i;
  }
  if(icb==0) {
    icbd = 0;
  } else {
    if(icbd<0) {
      if(l->myrank==0) {
	printf("not enough 2's in localGeom\n");
	PRINTV("physGeom:", "%i", l->physGeom, nd);
	PRINTV("rankGeom:", "%i", l->rankGeom, nd);
	PRINTV("localGeom:", "%i", l->localGeom, nd);
	PRINTV("outerGeom:", "%i", l->outerGeom, nd);
	PRINTV("innerGeom:", "%i", l->innerGeom, nd);
      }
      exit(-1);
    }
    icb = l->outerGeom[icbd]/2;
    if((icb&1)==0) {
      if(l->myrank==0) {
	printf("error in cb choice\n");
	PRINTV("physGeom:", "%i", l->physGeom, nd);
	PRINTV("rankGeom:", "%i", l->rankGeom, nd);
	PRINTV("localGeom:", "%i", l->localGeom, nd);
	PRINTV("outerGeom:", "%i", l->outerGeom, nd);
	PRINTV("innerGeom:", "%i", l->innerGeom, nd);
	printf("innerCb: %i\n", icb);
	printf("innerCbDir: %i\n", icbd);
      }
      exit(-1);
    }
  }
  l->physVol = pvol;
  l->nSites = lvol;
  l->nOdd = lvol/2;
  l->nEven = lvol - l->nOdd;
  l->nSitesOuter = ovol;
  l->nOddOuter = ovol/2;
  l->nEvenOuter = ovol - l->nOddOuter;
  l->nSitesInner = l->nSites/l->nSitesOuter;
  l->innerCb = icb;
  l->innerCbDir = icbd;
  if(l->myrank==0) {
    printf("#innerCb: %i\n", icb);
    printf("#innerCbDir: %i\n", icbd);
  }
}

static void
lex_x(int *x, int l, int *s, int ndim)
{
  for(int i=0; i<ndim; i++) {
  //for(int i=ndim-1; i>=0; --i) {
    x[i] = l % s[i];
    l = l / s[i];
  }
}

// x[0] is fastest
static int
lex_i(int *x, int *s, int *d, int ndim)
{
  int l = 0;
  for(int i=ndim-1; i>=0; --i) {
    int xx = x[i];
    if(d) xx /= d[i];
    l = l*s[i] + (xx%s[i]);
  }
  return l;
}

#if 0
// x[0] is slowest
static int
lexr_i(int *x, int *s, int *d, int ndim)
{
  int l = 0;
  for(int i=0; i<ndim; i++) {
    int xx = x[i];
    if(d) xx /= d[i];
    l = l*s[i] + (xx%s[i]);
  }
  return l;
}
#endif

void
layoutIndex(Layout *l, LayoutIndex *li, int coords[])
{
  int nd = l->nDim;
  int ri = lex_i(coords, l->rankGeom, l->localGeom, nd);
  int ii = lex_i(coords, l->innerGeom, l->outerGeom, nd);
  int ib = 0;
  for(int i=0; i<nd; i++) {
    int xi = coords[i]/l->outerGeom[i];
    int li = xi % l->innerGeom[i];
    ib += li * l->outerGeom[i];
  }
  ib &= 1;
  coords[l->innerCbDir] += l->innerCb * ib;
  int oi = lex_i(coords, l->outerGeom, NULL, nd);
  coords[l->innerCbDir] -= l->innerCb * ib;
  int p = 0;
  for(int i=0; i<nd; i++) p += coords[i];
  int oi2 = oi/2;
  if(p&1) oi2 = (oi+l->nSitesOuter)/2;
  li->rank = ri;
  li->index = oi2*l->nSitesInner + ii;
}

void
layoutCoord(Layout *l, int *coords, LayoutIndex *li)
{
  int nd = l->nDim;
  int cr[nd];
  lex_x(cr, li->rank, l->rankGeom, nd);
  int p = 0;
  int ll = li->index % l->nSitesInner;
  int ib = 0;
  for(int i=0; i<nd; i++) {
    int w = l->innerGeom[i];
    int wl = l->outerGeom[i];
    int k = ll % w;
    int c = l->localGeom[i]*cr[i] + k*wl;
    cr[i] = c;
    //printf("cr[%i]: %i\n", i, c);
    p += c;
    ll = ll / w;
    ib += k*wl;
  }
  ib &= 1;
  int ii = li->index / l->nSitesInner;
  if(ii>=l->nEvenOuter) {
    ii -= l->nEvenOuter;
    p++;
  }
  ii *= 2;
  for(int i=0; i<nd; i++) {
    int wl = l->outerGeom[i];
    int k = ii % wl;
    if(i==l->innerCbDir) k = (k + l->innerCb * ib)%wl;
    coords[i] = k;
    //printf("coords[%i]: %i\n", i, k);
    p += k;
    ii = ii / wl;
  }
  if(p&1) {
    for(int i=0; i<nd; i++) {
      int wl = l->outerGeom[i];
      if(i==l->innerCbDir) coords[i] = (coords[i] + l->innerCb * ib)%wl;
      coords[i]++;
      if(coords[i]>=wl) {
	coords[i] = 0;
	if(i==l->innerCbDir) coords[i] = (coords[i] + l->innerCb * ib)%wl;
      } else {
	if(i==l->innerCbDir) coords[i] = (coords[i] + l->innerCb * ib)%wl;
	break;
      }
    }
  }
  for(int i=0; i<nd; i++) coords[i] += cr[i];
  {
    LayoutIndex li2;
    layoutIndex(l, &li2, coords);
    if(li->rank!=li2.rank ||li->index!=li2.index) {
      printf("error: bad coord:\n");
      printf(" %i,%i -> %i %i %i %i -> %i,%i\n", li->rank, li->index,
	     coords[0],coords[1],coords[2],coords[3], li2.rank, li2.index);
      exit(-1);
    }
  }
}

void
layoutShift(Layout *l, LayoutIndex *li, LayoutIndex *li2, int *disp)
{
  int nd = l->nDim;
  int x[nd];
  layoutCoord(l, x, li2);
  for(int i=0; i<nd; i++) {
    x[i] = (x[i] + disp[i] + l->physGeom[i])%l->physGeom[i];
  }
  layoutIndex(l, li, x);
}
