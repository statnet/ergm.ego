/*  File src/changestats.c in package ergm.ego, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2015-2024 Statnet Commons
 */
#include "ergm_constants.h"
#include "ergm_changestat.h"

/*****************
 changestat: c_netsize_adj
*****************/
C_CHANGESTAT_FN(c_netsize_adj) {
  double edges_coef = INPUT_PARAM[0],
    mutual_coef = INPUT_PARAM[1],
    trties_coef = INPUT_PARAM[2],
    cyties_coef = INPUT_PARAM[3];

  if(edges_coef!=0)
    CHANGE_STAT[0] += edgestate ? - edges_coef : edges_coef;

  if(mutual_coef!=0)
    if(IS_OUTEDGE(head,tail)) /* otherwise, no change occurs */
      CHANGE_STAT[0] += edgestate ? - mutual_coef : mutual_coef;

  if(trties_coef!=0){
    int  echange, ochange;
    int L2th, L2tu, L2uh;
    double cumchange;
  
    cumchange=0.0;
    L2th=0;
    ochange = edgestate ? -1 : 0;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    EXEC_THROUGH_OUTEDGES(head,  e,  u, {
	if (GETWT(tail, u)){
	  L2tu=ochange;
	  /* step through inedges of u */
	  EXEC_THROUGH_INEDGES(u,  f,  v, {
	      if(GETWT(tail, v)){
		L2tu++;
		if(L2tu>0) {break;}
	      }
	    });
	  cumchange += (L2tu==0);
	}
      });
    /* step through inedges of head */
    
    EXEC_THROUGH_INEDGES(head,  e,  u, {
	if (GETWT(tail, u)){
	  L2th++;
	}
	if (GETWT(u, tail)){
	  L2uh=ochange;
	  /* step through outedges of u */
	  EXEC_THROUGH_OUTEDGES(u,  f,  v, {
	      if(GETWT(v, head)){
		L2uh++;
		if(L2uh>0) {break;}
	      }
	    });
	  cumchange += (L2uh==0) ;
	}
      });
  
    cumchange += (L2th>0) ;
    cumchange  = echange*cumchange;
    
    CHANGE_STAT[0] += cumchange*trties_coef;
  }

  if(cyties_coef!=0){
    int  echange, ochange;
    int L2th, L2tu, L2uh;
    double cumchange;
    
    /* *** don't forget tail -> head */    
    cumchange=0.0;
    L2th=0;
    ochange = edgestate ? -1 : 0;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    EXEC_THROUGH_OUTEDGES(head,  e,  u, {
      if (GETWT(u, tail)){
	L2tu=ochange;
	/* step through inedges of u */
	EXEC_THROUGH_INEDGES(u,  f,  v, {
	  if(GETWT(tail, v)){
	    L2tu++;
	    if(L2tu>0) {break;}
	  }
	});
	cumchange += (L2tu==0);
      }
    });
    /* step through outedges of head */
    
    EXEC_THROUGH_OUTEDGES(head,  e,  u, {
      if (GETWT(u, tail)){
	L2th++;
      }
      if (GETWT(u, tail)){
	L2uh=ochange;
	/* step through outedges of u */
	EXEC_THROUGH_OUTEDGES(u,  f,  v, {
	  if(GETWT(v, head)){
	    L2uh++;
	    if(L2uh>0) {break;}
	  }
	});
	cumchange += (L2uh==0) ;
      }
    });
    
    cumchange += (L2th>0) ;
    cumchange  = echange*cumchange;
    CHANGE_STAT[0] += cumchange*cyties_coef;
  }
}
