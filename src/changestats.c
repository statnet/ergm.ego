/*  File src/changestats.c in package ergm.ego, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2015-2022 Statnet Commons
 */
#include "ergm_constants.h"
#include "ergm_changestat.h"

#if ERGM_API_MAJOR >= 4
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
#else
/*****************
 changestat: d_netsize_adj
*****************/
D_CHANGESTAT_FN(d_netsize_adj) {
  double edges_coef = INPUT_PARAM[0],
    mutual_coef = INPUT_PARAM[1],
    trties_coef = INPUT_PARAM[2],
    cyties_coef = INPUT_PARAM[3];

  int i;
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    int  echange, ochange;
    Vertex tail, head, u, v;
    Edge e, f;
    Rboolean edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));

    if(edges_coef!=0)
      CHANGE_STAT[0] += edgestate ? - edges_coef : edges_coef;

    if(DIRECTED && mutual_coef!=0)
      if(IS_OUTEDGE(head,tail)) /* otherwise, no change occurs */
        CHANGE_STAT[0] += edgestate ? - mutual_coef : mutual_coef;

    ochange = edgestate ? -1 : 0;
    echange = 2*ochange + 1;

    if(DIRECTED){
      if(trties_coef!=0){
        int L2th, L2tu, L2uh;
        double cumchange;

        cumchange=0.0;
        L2th=0;
        /* step through outedges of head  */
        STEP_THROUGH_OUTEDGES(head,  e,  u){
          if (IS_OUTEDGE(tail, u)){
            L2tu=ochange;
            /* step through inedges of u */
            STEP_THROUGH_INEDGES(u,  f,  v){
              if(IS_OUTEDGE(tail, v)){
                L2tu++;
                if(L2tu>0) {break;}
              }
            }
            cumchange += (L2tu==0);
          }
        }
        /* step through inedges of head */

        STEP_THROUGH_INEDGES(head,  e,  u){
          if (IS_OUTEDGE(tail, u)){
            L2th++;
          }
          if (IS_OUTEDGE(u, tail)){
            L2uh=ochange;
            /* step through outedges of u */
            STEP_THROUGH_OUTEDGES(u,  f,  v){
              if(IS_OUTEDGE(v, head)){
                L2uh++;
                if(L2uh>0) {break;}
              }
            }
            cumchange += (L2uh==0) ;
          }
        }

        cumchange += (L2th>0) ;
        cumchange  = echange*cumchange;

        CHANGE_STAT[0] += cumchange*trties_coef;
      }

      if(cyties_coef!=0){
        int L2th, L2tu, L2uh;
        double cumchange;

        /* *** don't forget tail -> head */
        cumchange=0.0;
        L2th=0;
        /* step through outedges of head  */
        STEP_THROUGH_OUTEDGES(head,  e,  u){
          if (IS_OUTEDGE(u, tail)){
            L2tu=ochange;
            /* step through inedges of u */
            STEP_THROUGH_INEDGES(u,  f,  v){
              if(IS_OUTEDGE(tail, v)){
                L2tu++;
                if(L2tu>0) {break;}
              }
            }
            cumchange += (L2tu==0);
          }
        }
        /* step through outedges of head */

        STEP_THROUGH_OUTEDGES(head,  e,  u){
          if (IS_OUTEDGE(u, tail)){
            L2th++;
          }
          if (IS_OUTEDGE(u, tail)){
            L2uh=ochange;
            /* step through outedges of u */
            STEP_THROUGH_OUTEDGES(u,  f,  v){
              if(IS_OUTEDGE(v, head)){
                L2uh++;
                if(L2uh>0) {break;}
              }
            }
            cumchange += (L2uh==0) ;
          }
        }

        cumchange += (L2th>0) ;
        cumchange  = echange*cumchange;
        CHANGE_STAT[0] += cumchange*cyties_coef;
      }
    }else{
      if(trties_coef!=0){
        Edge e, f;
        int L2th, L2tu, L2uh;
        Vertex u, v;
        double cumchange;

        /* *** don't forget tail -> head */
        cumchange=0.0;
        L2th=0;
        /* step through outedges of head  */
        STEP_THROUGH_OUTEDGES(head, e, u){
          if (IS_UNDIRECTED_EDGE(u, tail)){
            L2th++;
            L2tu=ochange;
            L2uh=ochange;
            /* step through outedges of u */
            STEP_THROUGH_OUTEDGES(u, f, v){
              if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
              if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
            }
            /* step through inedges of u */
            STEP_THROUGH_INEDGES(u, f, v){
              if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
              if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
            }
            cumchange += (L2tu==0) + (L2uh==0);
          }
        }
        /* step through inedges of head */
        STEP_THROUGH_INEDGES(head, e, u){
          if (IS_UNDIRECTED_EDGE(u, tail)){
            L2th++;
            L2tu=ochange;
            L2uh=ochange;
            /* step through outedges of u */
            STEP_THROUGH_OUTEDGES(u, f, v){
              if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
              if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
            }
            /* step through inedges of u */
            STEP_THROUGH_INEDGES(u, f, v){
              if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
              if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
            }
            cumchange += (L2tu==0) + (L2uh==0);
          }
        }

        cumchange += (L2th!=0);
        cumchange *= echange;
        CHANGE_STAT[0] += cumchange*trties_coef;
      }
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}
#endif
