#include "VecMat.h"

//NXM�ľ���A�е���row0��row1Ԫ�ػ���
template<class T> inline void swaprows(T* A, int row0, int row1, int N, int M)
{
    T temp=0;
	int s0=row0*M,s1=row1*M;
	for(int i=0;i<M;i++)
	{
		temp=A[s0];
		A[s0++]=A[s1];
		A[s1++]=temp;
	}
}
//	gjelim 
//A,B:NXN,WA:NX2N
void M_Inverse(double* B, const double* A, int N, double* WA)
{
	int row,col,dindex, N2=2*N;
	int swa=0,sa=0;
    for (row=0; row<N; ++row)
	{
        for (col=0; col<N; ++col)
            WA[swa++]=A[sa++];
        for (col=N; col<N2; ++col)
		{
			if(row==col-N) WA[swa++]=1.0;
			else WA[swa++]=0.0;
		}
    }
    //	perform forward elimination to get WA in row-echelon form
    for (dindex=0; dindex<N; ++dindex)
	{
		int ss=dindex*N2+dindex;
        //	run along diagonal, swapping rows to move zeros in working position 
        //	(along the diagonal) downwards
        if ( (dindex==(N-1)) && (WA[ss]==0)) return; //  no solution
	//	else if (WA[ss]==0) swaprows(WA, dindex, dindex+1,N,N2);
		int maxrow=dindex;
		double maxvalue=fabs(WA[ss]);
		for(int kk=dindex+1;kk<N;kk++)
		{
			if (fabs(WA[kk*N2+dindex])>maxvalue) {maxvalue=fabs(WA[kk*N2+dindex]);maxrow=kk;}
		}
		if (maxrow!=dindex) swaprows(WA, dindex, maxrow,N,N2);
        //	divide working row by value of working position to get a 1 on the
        //	diagonal
        if (WA[ss] == 0.0) return;
		else
		{
            double tempval=WA[ss];
			int index=dindex*N2;
            for (col=0; col<N2; ++col) WA[index++]/=tempval;
        }

        //	eliminate value below working position by subtracting a multiple of 
        //	the current row
        for (row=dindex+1; row<N; ++row)
		{
            double wval=WA[row*N2+dindex];
            for (col=0; col<N2; ++col) WA[row*N2+col]-=wval*WA[dindex*N2+col];
        }
    }

    //	backward substitution steps
    for (dindex=N-1; dindex>=0; --dindex)
	{
        //	eliminate value above working position by subtracting a multiple of 
        //	the current row
        for (row=dindex-1; row>=0; --row)
		{
            double wval=WA[row*N2+dindex];
            for (col=0; col<N2; ++col)
                WA[row*N2+col]-=wval*WA[dindex*N2+col];
        }
    }

    //	assign result to replace B
    for (row=0; row<N; ++row)
	{
        for (col=0; col<N; ++col)
			B[row*N+col]=WA[row*N2+col+N];
    }
}
/*
void swaprows(double** arr, long row0, long row1)
 {
    double* temp;
    temp=arr[row0];
    arr[row0]=arr[row1];
    arr[row1]=temp;
}

void M_Inverse(double** B, double** A, int dim, double** WA)
{
	int row,col,dindex;
    //	augment A array with B array and store in WA
//    double** WA=new double*[dim];
 //   for (int row=0; row<dim; ++row)
 //       WA[row]=new double[dim+dim];
    for (row=0; row<dim; ++row)
	{
        for (col=0; col<dim; ++col)
            WA[row][col]=A[row][col];
        for (col=dim; col<dim+dim; ++col)
		{
			if(row==col-dim) WA[row][col]=1.0;
			else WA[row][col]=0.0;
		}
    }

    //	perform forward elimination to get WA in row-echelon form
    for (dindex=0; dindex<dim; ++dindex)
	{
        //	run along diagonal, swapping rows to move zeros in working position 
        //	(along the diagonal) downwards
        if ( (dindex==(dim-1)) && (WA[dindex][dindex]==0)) return; //  no solution
		int maxrow=dindex;
		double maxvalue=fabs(WA[dindex][dindex]);
		for(int kk=dindex+1;kk<dim;kk++)
		{
			if (fabs(WA[kk][dindex])>maxvalue) {maxvalue=fabs(WA[kk][dindex]);maxrow=kk;}
		}
		if (maxrow!=dindex) swaprows(WA, dindex, maxrow);
        //	divide working row by value of working position to get a 1 on the
        //	diagonal
        if (WA[dindex][dindex] == 0.0) return;
		else
		{
            double tempval=WA[dindex][dindex];
            for (col=0; col<dim+dim; ++col) WA[dindex][col]/=tempval;
        }

        //	eliminate value below working position by subtracting a multiple of 
        //	the current row
        for (row=dindex+1; row<dim; ++row)
		{
            double wval=WA[row][dindex];
            for (col=0; col<dim+dim; ++col) WA[row][col]-=wval*WA[dindex][col];
        }
    }

    //	backward substitution steps
    for (dindex=dim-1; dindex>=0; --dindex)
	{
        //	eliminate value above working position by subtracting a multiple of 
        //	the current row
        for (row=dindex-1; row>=0; --row)
		{
            double wval=WA[row][dindex];
            for (col=0; col<dim+dim; ++col)
                WA[row][col]-=wval*WA[dindex][col];
        }
    }

    //	assign result to replace B
    for (row=0; row<dim; ++row)
	{
        for (col=0; col<dim; ++col)
			B[row][col]=WA[row][col+dim];
    }

//    for (int row=0; row<dim; ++row)
//       delete[] WA[row];
//    delete[] WA;
}*/