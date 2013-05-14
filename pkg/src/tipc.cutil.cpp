#include <iostream>
#include <fstream>
#include <cstdlib>
#include "tipc.cutil.h"

static inline int cmpi(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}


SEXP tipc_tabulate_after_sample(SEXP tpc, SEXP sid)
{
	int i=0, j= 0, m=0, n= length(tpc), nsid= length(sid), sampledtree= 0;
	int tipc_m_ncol= 4, tipc_m_nrow= 0;//assume four status levels
	int *xid= NULL, *yid=NULL, *ystatus= NULL, *xsid= INTEGER(sid);
	int *tipc_tab= NULL, *xtipc_m= NULL, *ytipc_m= NULL;
	SEXP tpc_trees= VECTOR_ELT(tpc,0);
	SEXP tree= NULL, tree_nodes= NULL, tree_nodes_id= NULL, tree_nodes_status= NULL;
	SEXP tipc_m= NULL;

	n= length(tpc_trees);
	xtipc_m= tipc_tab= NEW_ZERO_ARY(int, tipc_m_ncol*n);//create large enough array


	for(i= 0; i<n; i++)//loop over all tpc_trees
	{
		//get pointers to parts of data frame
		tree= VECTOR_ELT(tpc_trees,i);
		tree_nodes= VECTOR_ELT(tree,0);
		tree_nodes_id= VECTOR_ELT(tree_nodes,0);
		tree_nodes_status= VECTOR_ELT(tree_nodes,1);

		m= length(tree_nodes_id);
		xid= INTEGER(tree_nodes_id);
		ystatus= INTEGER(tree_nodes_status);
		for(j=0, sampledtree= 0; 		j<m; 		j++, xid++, ystatus++)//for all nodes in node data frame
		{
			yid = (int*)bsearch(xid, xsid, nsid, sizeof(int), cmpi);
			if(!j && !yid)	break;//first node is root and must be sampled; otherwise abort
			if(yid)
			{
				(*(xtipc_m+*ystatus-1))++;
				sampledtree++;
			}
			/*if(yid)
				printf ("%d is in the array. status is %d\n",*yid, *ystatus);
			else
				printf ("%d is not in the array. status is %d\n",*xid, *ystatus);
			*/
		}
		if(sampledtree)
		{
			xtipc_m+=tipc_m_ncol;
			tipc_m_nrow++;
		}
	}
	PROTECT(tipc_m = allocMatrix(INTSXP, tipc_m_ncol, tipc_m_nrow ));//increments first over rows
	for(i= tipc_m_ncol*tipc_m_nrow, xtipc_m= INTEGER(tipc_m), ytipc_m= tipc_tab; 	i--; 	*xtipc_m++=*ytipc_m++ );

	DELETE(tipc_tab);
	UNPROTECT(1);
	return tipc_m;
}

