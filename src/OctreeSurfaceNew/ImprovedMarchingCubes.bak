#define BIG 1
#define SMALL 2
#define DEFAULT_ERR		0.0001
//#define DEFAULT_IVAL	0.0001
#define DEFAULT_IVAL	0.0001
//#define DEFAULT_IVAL	-5.5001
#define SIGN_POSITIVE 1
#define SIGN_NEGATIVE -1

class Scaffold
{
public :
	int nvert;
	int ntet;
	int npyramid;
	int ncube;

	double vertices[17][3];
	double func[17];

	int tet[24][4];
	int pyramid[6][5];


	Scaffold()
	{
		nvert=0;
		ntet=0;
		npyramid=0;
		ncube;
	}

	~Scaffold()
	{
	}

	void clear()
	{
		nvert=0;
		ntet=0;
		npyramid=0;
	}
	void addVert(double pos[3])
	{
		int i;
		for (i=0;i<3;i++) {
			vertices[nvert][i]=pos[i];
		}
		nvert++;
	}

	void addTet(int input_tet[4])
	{
		int i;
		for (i=0;i<4;i++) tet[ntet][i]=input_tet[i];
		ntet++;
	}

	void addPyramid(int input_pyramid[5])
	{
		int i;
		for (i=0;i<5;i++) pyramid[npyramid][i]=input_pyramid[i];
		npyramid++;
	}

};

class CellQueue {
   public:
      // constructor/destructor
      inline CellQueue(int size=100);
      ~CellQueue();

      // add item to the queue
      inline void Add(unsigned int cell);

      // remove and return the first item in queue
      inline int  Get(int &cell);

      // return the first item in queue
      inline int  Peek(int &cell);

      // remove the first item in queue
      inline void Pop();

      // reset to empty
      void Reset(void) { nel = 0; }

      // check if queue is empty
      int  Empty(void) { return(nel == 0); }

   protected:

   private:
      int nel;
      int cellsize;  /* # of elements in cell array   */
      int start;
      unsigned int *cells;
};

//------------------------------------------------------------------------
//
// CellQueue() - create a new cell queue with elements of specified size
//
//------------------------------------------------------------------------
inline CellQueue::CellQueue(int size)
{
   nel    = 0;
   start  = 0;
   cellsize = size;
   cells    = (unsigned int *)malloc(sizeof(unsigned int) * cellsize);
}


//------------------------------------------------------------------------
//
// ~CellQueue() - free storage
//
//------------------------------------------------------------------------
inline CellQueue::~CellQueue()
{
   if (cells != NULL)
      free(cells);
}


//------------------------------------------------------------------------
//
// Add() - add an item to the queue
//
//------------------------------------------------------------------------
inline void
CellQueue::Add(unsigned int c)
{
   int n;
   int oldsize;
   int atend;

   n = nel++;

   // resize the queue if needed
   if (nel > cellsize) {
      oldsize = cellsize;
      cellsize *= 2;
      cells = (unsigned int *)realloc(cells, sizeof(int) * cellsize);

      // move everything from 'start' to the end
      if (start != 0) {
         atend = oldsize - start;
         memmove(&cells[cellsize-atend], &cells[start], sizeof(unsigned int)*atend);
         start = cellsize-atend;
      }
   }

   n += start;
   if (n >= cellsize)
      n-=cellsize;

   cells[n] = c;
}

//------------------------------------------------------------------------
//
// Get() - return the top item from the queue
//
//------------------------------------------------------------------------
inline int
CellQueue::Get(int &c)
{
   if (Peek(c) == -1)
      return(-1);
   Pop();
   return(1);
}

//------------------------------------------------------------------------
//
// Peek() - return the top item, but don't remove it
//
//------------------------------------------------------------------------
inline int
CellQueue::Peek(int &c)
{
   if (nel == 0)
      return(-1);

   c = cells[start];

   return(1);
}

//------------------------------------------------------------------------
//
// Pop() - delete the top item in the queue
//
//------------------------------------------------------------------------
inline void
CellQueue::Pop(void)
{
   start++;
   if (start == cellsize)
      start=0;
   nel--;
}



int shared_face_edge[6][6][2];

int vtx_connectivity[8][3]={{1,2,4} , {0,3,5} , {0,3,6} , {1,2,7} , {0,5,6} , {1,4,7} , {2,4,7} , {3,5,6}};

int tet_edge[6][2] = {
	{0,1},
	{0,2},
	{0,3},
	{1,2},
	{2,3},
	{3,1}
};

int tet_triangle[16][7] = {
	{0,-1,-1,-1,-1,-1,-1} ,
	{1, 0, 1, 2,-1,-1,-1},
	{1, 0, 3, 5,-1,-1,-1},
	{2, 1, 3, 5, 1, 5, 2},
	{1, 1, 3, 4,-1,-1,-1},
	{2, 0, 3, 4, 4, 0, 2},
	{2, 1, 4, 5, 0, 1, 5},
	{1, 2, 4, 5,-1,-1,-1},
	{1, 2, 4, 5,-1,-1,-1},
	{2, 1, 4, 5, 0, 1, 5},
	{2, 0, 3, 4, 4, 0, 2},
	{1, 1, 3, 4,-1,-1,-1},
	{2, 1, 3, 5, 1, 5, 2},
	{1, 0, 3, 5,-1,-1,-1},
	{1, 0, 1, 2,-1,-1,-1},
	{0,-1,-1,-1,-1,-1,-1} 
};

int pyramid_edge[8][2] = {
	{0,1},
	{0,2},
	{0,3},
	{0,4},
	{1,2},
	{2,4},
	{4,3},
	{3,1}
};

int pyramid_triangle[32][10] = {
	{0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1, 0, 4, 7,-1,-1,-1,-1,-1,-1},
	{1, 1, 5, 4,-1,-1,-1,-1,-1,-1},
	{2, 0, 7, 5, 0, 5, 1,-1,-1,-1},
	{1, 6, 7, 2,-1,-1,-1,-1,-1,-1},
	{2, 0, 2, 4, 2, 4, 6,-1,-1,-1},
	{0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3, 5, 6, 2, 5, 2, 0, 5, 0, 1},
	{1, 5, 6, 3,-1,-1,-1,-1,-1,-1},
	{0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2, 1, 4, 6, 1, 6, 3,-1,-1,-1},
	{3, 7, 6, 3, 7, 3, 1, 7, 1, 0},
	{2, 7, 2, 5, 2, 5, 3,-1,-1,-1},
	{3, 0, 2, 4, 4, 2, 3, 4, 3, 5},
	{3, 4, 7, 2, 4, 2, 3, 4, 3, 1},
	{2, 0, 2, 3, 0, 3, 1,-1,-1,-1},
	
	{2, 0, 2, 3, 0, 3, 1,-1,-1,-1},
	{3, 4, 7, 2, 4, 2, 3, 4, 3, 1},
	{3, 0, 2, 4, 4, 2, 3, 4, 3, 5},
	{2, 7, 2, 5, 2, 5, 3,-1,-1,-1},
	{3, 7, 6, 3, 7, 3, 1, 7, 1, 0},
	{2, 1, 4, 6, 1, 6, 3,-1,-1,-1},
	{0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1, 5, 6, 3,-1,-1,-1,-1,-1,-1},
	{3, 5, 6, 2, 5, 2, 0, 5, 0, 1},
	{0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2, 0, 2, 4, 2, 4, 6,-1,-1,-1},
	{1, 6, 7, 2,-1,-1,-1,-1,-1,-1},
	{2, 0, 7, 5, 0, 5, 1,-1,-1,-1},
	{1, 1, 5, 4,-1,-1,-1,-1,-1,-1},
	{1, 0, 4, 7,-1,-1,-1,-1,-1,-1},
	{0,-1,-1,-1,-1,-1,-1,-1,-1,-1}
};
	


double v_pos[8][3] = {
	{0,0,0},
	{1,0,0},
	{0,1,0},
	{1,1,0},
	{0,0,1},
	{1,0,1},
	{0,1,1},
	{1,1,1}};

double center_pos[3] = { 0.5f,0.5f,0.5f };
#define CUBE_VTX_NO 8
#define FACE_NO 6
#define BODY_NO 2
#define BODY_CENTER 16
#define BODY_SADDLE1 14
#define BODY_SADDLE2 15

int face_set[6][4] = {
	{0,1,2,3},
	{4,5,0,1},
	{1,5,3,7},
	{2,3,6,7},
	{0,4,2,6},
	{4,5,6,7}
};

int sixfold_tet[6][4] = {
	{0,1,3,2},
	{4,5,1,0},
	{5,7,3,1},
	{6,7,3,2},
	{4,6,2,0},
	{4,5,7,6}
};

int face_line[4][2] = { {0,1}, {1,3}, {3,2}, {2,0} };


int Vtx_Edge(int vid, int fi, int k)
{
  int i;
  for (i=0;i<4;i++) {
    if (k==0 && face_set[fi][face_line[i][0]]==vid)  return face_set[fi][face_line[i][1]];
    if (k==1 && face_set[fi][face_line[i][1]]==vid)  return face_set[fi][face_line[i][0]];
  }
  return -1;
  //assert(0);
}

int face_diagonal(int fi, int vid)
{
  int i,v0,v1		;
  for (i=0;i<4;i++) {
    if (face_set[fi][face_line[i][0]]==vid) v0=face_set[fi][face_line[i][1]];
    if (face_set[fi][face_line[i][1]]==vid) v1=face_set[fi][face_line[i][0]];
  }
  
  for (i=0;i<4;i++) {
    if (face_set[fi][i]!=vid && face_set[fi][i]!=v0 && face_set[fi][i]!=v1) return face_set[fi][i];
  }
  return -1;
  //assert(0);
}


int get_shared_vtx(int f1, int f2, int f3) 
{
  int v,i,j,k;
  for (v=0;v<8;v++) {
    for (i=0;i<4;i++) {
      for (j=0;j<4;j++) {
	for (k=0;k<4;k++) {
	  if (v==face_set[f1][i] && v==face_set[f2][j] && v==face_set[f3][k])  return v;
	}
      }
    }
  }  
  return -1;
}

int getPositiveNo(double* val, double isovalue)
{
  int i;
  int count=0;
  
  for (i=0;i<8;i++) {
    if (val[i]>=isovalue) count++;
  }
  return count;
}

int vertex_sign(double a,double isovalue)
{
	if (a>=isovalue) return SIGN_POSITIVE;
	else return SIGN_NEGATIVE;
}


bool is_contain(int fi , int vid)
{
  int i;
  for (i=0;i<4;i++) {
    if (face_set[fi][i]==vid) return true;
  } 
  return false;
}
void get_vtx_adj_face(int vid, int* face_array)
{
  int i,count;
  
  count=0;
  
  for (i=0;i<6;i++) {
    if (is_contain(i,vid)) face_array[count++]=i;
  }
  //assert(count==3);
}

double get_distance(double* s1, int face_id)
{
  switch (face_id) {
  case 0 :
    return s1[2];
  case 1 :
    return s1[1];
  case 2:
    return 1-s1[0];
  case 3:
    return 1-s1[1];
  case 4 :
    return s1[0];
  case 5 :
    return 1-s1[2];
  }
  return 1;
}

bool is_common_edge(int f1, int f2)
{
  
  if (shared_face_edge[f1][f2][0]==-1 || 
      shared_face_edge[f1][f2][1]==-1) 
    return false;
  return true;
}


double getTriVal(double val[8], double x  , double y  , double z  , int res)
{
  double x_ratio,y_ratio,z_ratio;
  double temp1,temp2,temp3,temp4,temp5,temp6;
  
  x_ratio=((double)(res-x))/((double)res);
  y_ratio=((double)(res-y))/((double)res);
  z_ratio=((double)(res-z))/((double)res);
  
  temp1 = val[1] + (val[0]-val[1])*x_ratio;
  temp2 = val[3] + (val[2]-val[3])*x_ratio;
  temp3 = val[5] + (val[4]-val[5])*x_ratio;
  temp4 = val[7] + (val[6]-val[7])*x_ratio;
  
  temp5 = temp2  + (temp1-temp2)*y_ratio;
  temp6 = temp4  + (temp3-temp4)*y_ratio;
  
  return temp6  + (temp5-temp6)*z_ratio;
}

void getBigSmallVertex(Scaffold& scaf, double* val, double big_saddle_val,double small_saddle_val, int* big_vid, int* small_vid)
{	
  int i,j;
  int count;
  double temp_val;
  int temp_idx;
  
  int saddle_idx[6];
  double saddle_val[6];
  
  for (i=0;i<6;i++) {
    saddle_idx[i]=i;
    saddle_val[i]=scaf.func[CUBE_VTX_NO+i];
  }
  
  for (i=0;i<5;i++) {
    for (j=0;j<5-i;j++) {
      if (saddle_val[j]>=saddle_val[j+1]) {
	temp_val=saddle_val[j];
	saddle_val[j]=saddle_val[j+1];
	saddle_val[j+1]=temp_val;
	
	temp_idx=saddle_idx[j];
	saddle_idx[j]=saddle_idx[j+1];
	saddle_idx[j+1]=temp_idx;
      }
    }
  }
  
  *big_vid=get_shared_vtx(saddle_idx[0],saddle_idx[1],saddle_idx[2]);
  *small_vid=get_shared_vtx(saddle_idx[3],saddle_idx[4],saddle_idx[5]);
}

int getBiggestVertex(double* val)
{
  int i;
  int max_idx=0;
  double max_val=val[0];
  
  for (i=1;i<8;i++) {
    if (val[i]>max_val) {
      max_idx=i;
      max_val=val[i];
    }
  }
  return max_idx;
}

int getSmallestVertex(double* val)
{
  int i;
  int min_idx=0;
  double min_val=val[0];
  
  for (i=1;i<8;i++) {
    if (val[i]<min_val) {
      min_idx=i;
      min_val=val[i];
    }
  }
  return min_idx;
}

int getThirdBiggestFaceSaddle(double* func,int big_vid)
{
	int i,j;
	int count;
	double temp_val;
	int temp_idx;
	
	int saddle_idx[6];
	double saddle_val[6];

	for (i=0;i<6;i++) {
		saddle_idx[i]=i;
		saddle_val[i]=func[CUBE_VTX_NO+i];		
	}
	
	for (i=0;i<5;i++) {
		for (j=0;j<5-i;j++) {
			if (saddle_val[j]>=saddle_val[j+1]) {
				temp_val=saddle_val[j];
				saddle_val[j]=saddle_val[j+1];
				saddle_val[j+1]=temp_val;
				temp_idx=saddle_idx[j];
				saddle_idx[j]=saddle_idx[j+1];
				saddle_idx[j+1]=temp_idx;
			}
		}
	}
	return saddle_idx[3];
}


int getThirdSmallestFaceSaddle(double* func,int small_vid)
{
	int i,j;
	int count;
	double temp_val;
	int temp_idx;
	
	int saddle_idx[6];
	double saddle_val[6];

	for (i=0;i<6;i++) {
		saddle_idx[i]=i;
		saddle_val[i]=func[CUBE_VTX_NO+i];
	}
	
	for (i=0;i<5;i++) {
		for (j=0;j<5-i;j++) {
			if (saddle_val[j]>=saddle_val[j+1]) {
				temp_val=saddle_val[j];
				saddle_val[j]=saddle_val[j+1];
				saddle_val[j+1]=temp_val;
				
				temp_idx=saddle_idx[j];
				saddle_idx[j]=saddle_idx[j+1];
				saddle_idx[j+1]=temp_idx;
			}
		}
	}
	return saddle_idx[2];
}


void determine_two_face(double* val, double* body_saddle1, double* body_saddle2, int* face_id1 , int* face_id2, int big_vid, int small_vid, int no_bsaddle)
{
  int i;
  double d[6];
  
  int face_pair[3][2]={
    {0,5},
    {1,3},
    {2,4}};
  
  double min_distance=10;
  int min_idx;
  int big_fsaddle;
  
  if (no_bsaddle==0) {
    
    
    big_fsaddle=getThirdBiggestFaceSaddle(val,big_vid);
    *face_id1=big_fsaddle;
    min_idx=big_fsaddle;
    
    
    for (i=0;i<3;i++) {
      if (face_pair[i][0]==min_idx) {
	*face_id2=face_pair[i][1];
      } else if (face_pair[i][1]==min_idx) {
				*face_id2=face_pair[i][0];
      }
    }
    
    /*
      
      if (val[CUBE_VTX_NO]>val[CUBE_VTX_NO+5]) {
      *face_id1=5; *face_id2=0;
      } else {
      *face_id1=0; *face_id2=5;
      }
      return;
    */
  }
  
  int vtx_face_array[3];
  if (no_bsaddle==1) {
    if (val[BODY_SADDLE2]<val[BODY_SADDLE1]) {
      // big_saddle
      
      
      get_vtx_adj_face(small_vid,vtx_face_array);
      for (i=0;i<3;i++) {
	d[i]=get_distance(body_saddle1,vtx_face_array[i]);
	if (min_distance>d[i]) {
	  min_idx=vtx_face_array[i];
	  min_distance=d[i];
	}
      }
      
      *face_id1=min_idx;
      for (i=0;i<3;i++) {
	if (face_pair[i][0]==min_idx) {
	  *face_id2=face_pair[i][1];
	} else if (face_pair[i][1]==min_idx) {
	  *face_id2=face_pair[i][0];
	}
      }
      
      
    } else {
      // small_saddle
      
      get_vtx_adj_face(big_vid,vtx_face_array);
      for (i=0;i<3;i++) {
	d[i]=get_distance(body_saddle1,vtx_face_array[i]);
	if (min_distance>d[i]) {
	  min_idx=vtx_face_array[i];
	  min_distance=d[i];
	}
      }
      
      *face_id2=min_idx;
      for (i=0;i<3;i++) {
	if (face_pair[i][0]==min_idx) {
	  *face_id1=face_pair[i][1];
	} else if (face_pair[i][1]==min_idx) {
	  *face_id1=face_pair[i][0];
	}
      }
      
    }
    return;
  }
  
  for (i=0;i<3;i++) {
    d[2*i]=get_distance(body_saddle1,face_pair[i][0])+get_distance(body_saddle2,face_pair[i][1]);
    if (min_distance>d[2*i]) {
      min_idx=2*i;
      min_distance=d[2*i];
    }
    
    d[2*i+1]=get_distance(body_saddle1,face_pair[i][1])+get_distance(body_saddle2,face_pair[i][0]);		
    if (min_distance>d[2*i+1]) {
      min_idx=2*i+1;
      min_distance=d[2*i+1];
    }
  }
  
  switch (min_idx) {
  case 0 :
    *face_id1=0; *face_id2=5;
    break;
  case 1 :
    *face_id1=5; *face_id2=0;
    break;
    
  case 2:
    *face_id1=1; *face_id2=3;
    break;
  case 3 :
    *face_id1=3; *face_id2=1;
    break;
    
  case 4 :
    *face_id1=2; *face_id2=4;
    break;
  case 5 :
    *face_id1=4; *face_id2=2;
    break;
  }
}




bool isInCube(double* temp_vtx) 
{
  int i;
  for (i=0;i<3;i++) {
    if (temp_vtx[i]<=0 || temp_vtx[i]>=1) return false;
  }
  return true;
}

bool isFaceAmbiguous(double* val, double isovalue)
{
  int i;
  int code[4];
  int count=0;
  
  for (i=0;i<4;i++) {
    if (val[i]>=isovalue) {
      code[i]=1;
      count++;
    } else code[i]=0;
  }
  if (count==2 && code[0]==code[3]) return true;
  else return false;  
}


bool is_ambiguous_case(double* val,double iso_val,int Case)
{
  CellQueue vtx_queue;
  //int vtx_connectivity[8][3]={{1,3,4} , {0,2,5} , {1,3,6} , {0,2,7} , {0,5,7} , {1,4,6} , {2,5,7} , {3,4,6}};
  int vtx_connectivity[8][3]={{1,2,4} , {0,3,5} , {0,3,6} , {1,2,7} , {0,5,6} , {1,4,7} , {2,4,7} , {3,5,6}};
  int vtx_idx=-1;
  
  int code[8] = {0,0,0,0,0,0,0,0} ;
  int i;
  for (i = 0 ; i < 8 ; i++) {
    if (Case==0) {
      if (val[i] >= iso_val) {
	code[i] = 1;
	vtx_idx=i;
      }
    } else {
      if (val[i] < iso_val) {
	code[i] = 1;
	vtx_idx=i;
      }
    }
  }
  
  if (vtx_idx==-1) return false;
  else {
    vtx_queue.Add(vtx_idx);
    code[vtx_idx]=0;
    
    while (vtx_queue.Get(vtx_idx)>=0) {
      for (i = 0 ; i < 3 ; i++) {
	if (code[vtx_connectivity[vtx_idx][i]]==1) {
	  vtx_queue.Add(vtx_connectivity[vtx_idx][i]);
	  code[vtx_connectivity[vtx_idx][i]]=0;
	}
      }
    }
  }
  
  for (i = 0 ; i < 8 ; i++) {
    if (code[i]==1) return true;
  }
  return false;
}


// this function decides whether an isosurface in the cell specified by
// its value is ambiguous or not.
// if a cell has an ambiguous cell, you must refine the cell until
// isosurfaces in the cells are not ambiguous
bool isType1(double* val,double iso_val)
{
        bool flag_val_case1,flag_val_case2;
        flag_val_case1=is_ambiguous_case(val,iso_val,0);
        flag_val_case2=is_ambiguous_case(val,iso_val,1);
        if (flag_val_case1==true || flag_val_case2==true) return false;
        else return true;
}


bool isType2(double* val,double isovalue, int bno)
{
	int i,j;
	double face_val[4];
	int ambiguous_face_no=0;

	for (i=0;i<6;i++) {
		for (j=0;j<4;j++) {
			face_val[j]=val[face_set[i][j]];
		}
		if (isFaceAmbiguous(face_val,isovalue)) {
			ambiguous_face_no++;
		}
	}
	if (ambiguous_face_no<3 && ambiguous_face_no>=0) return true;
	//else if (ambiguous_face_no==0 && bno==1) return true;
	else return false;
}

bool isType3(double* val,double isovalue)
{
	int i,j;
	double face_val[4];
	int ambiguous_face_no=0;

	for (i=0;i<6;i++) {
		for (j=0;j<4;j++) {
			face_val[j]=val[face_set[i][j]];
		}
		if (isFaceAmbiguous(face_val,isovalue)) {
			ambiguous_face_no++;
		}
	}
	if (ambiguous_face_no==3) return true;
	else return false;
}



bool computeFaceSaddle(double* fval, double* fsaddle,int face_idx)
{
  double p,q,r,s,a,b;
  
  
  p=fval[0];
  q=-fval[0]+fval[1];
  r=-fval[0]+fval[2];
  s=fval[0]-fval[1]-fval[2]+fval[3];
  
  a=-q/s;
  b=-r/s;
  
  switch (face_idx) {
  case 0 :
    fsaddle[0]=b;
    fsaddle[1]=a;
    fsaddle[2]=0.0f;
    break;
  case 1:
    fsaddle[0]=b;
    fsaddle[1]=0.0f;
    fsaddle[2]=1-a;
    break;
  case 2:
    fsaddle[0]=1;
    fsaddle[1]=a;
    fsaddle[2]=b;
    break;
    
  case 3 :
    fsaddle[0]=b;
    fsaddle[1]=1.0f;
    fsaddle[2]=a;
    break;
    
  case 4:
    fsaddle[0]=0.0f;
    fsaddle[1]=a;
    fsaddle[2]=b;
    break;
    
  case 5:
    fsaddle[0]=b;
    fsaddle[1]=a;
    fsaddle[2]=1.0f;
    break;
  }
  
  if ((a>=0 && a<=1) && (b>=0 && b<=1)) return true;
  else return false; 
  
}


int computeBodySaddle(double* bval, double* bsaddle1, double* bsaddle2)
//int computeBodySaddle(float* bval,float iso_val,vector<char> &bsaddlesigns)
{
	double a,b,c,d,e,f,g,h,k0,v1,v2,z1,z2;
	double aebd, afbc, cdag;
	//float bsaddle1[3],bsaddle2[3];
	
	a = bval[1]+bval[2]+bval[3]+bval[7]-bval[6]-bval[4]-bval[5]-bval[0];
	b = bval[4]-bval[2]-bval[1]+bval[0];
	c = bval[6] - bval[3] - bval[1] + bval[0];
	d = bval[5] -bval[3] - bval[2] + bval[0];
	e = bval[1] - bval[0];
	f = bval[2] - bval[0];
	g = bval[3] - bval[0];
	h = bval[0];
	
	aebd = a*e-b*d; 
	afbc = a*f-b*c;
	cdag = c*d-a*g;
	
	if(fabs(aebd)<10e-6)
	{
		cout<<"\n Warning: aebd may be zero. Body saddle might run into problems."<<endl;
	}
//	cout<<"aebd: "<<aebd<<endl;
	k0 = afbc*cdag/aebd;
	if(k0<0.0)
	{
      cout<<"\n Warning: Imaginary body saddle found.  "<<"k0= "<<k0<<endl;
	//	bsaddlesigns.push_back(1); // bsaddle sign assumed to be 0.
	//	return 1; // One body saddle by default.
      //return -1; 
	    return 0; 
    }
	else 
	{
		z1 = (-1*c + sqrt(k0))/a;
		z2 = (-1*c - sqrt(k0))/a;
		
		bsaddle1[0] = z1;
		bsaddle1[1] = -1*(d*z1+g)/(a*z1+c);
		bsaddle1[2] = -1*(b*z1+f)/(a*z1+c);
		v1 = e*z1 + b - (d*z1+g)*(b*z1+f)/(a*z1+c);
		
		bsaddle2[0] = z2;
		bsaddle2[1] = -1*(d*z2+g)/(a*z2+c);
		bsaddle2[2] = -1*(b*z2+f)/(a*z2+c);
		v2 = e*z2 + b - (d*z2+g)*(b*z2+f)/(a*z2+c);
		
		float temp_vtx[3],temp_val;
		
		if (isInCube(bsaddle1)) {
		    if (isInCube(bsaddle2)) {
		      if (v1 < v2) {
		        for (int i=0;i<3;i++) temp_vtx[i]=bsaddle2[i];
		        for (int i=0;i<3;i++) bsaddle2[i]=bsaddle1[i];
		        for (int i=0;i<3;i++) bsaddle1[i]=temp_vtx[i];
				temp_val=v1;
				v1 = v2;
				v2 = temp_val;
		      }
		      //if(v1>=iso_val) bsaddlesigns.push_back(0); else bsaddlesigns.push_back(1);
		      //if(v2>=iso_val) bsaddlesigns.push_back(0); else bsaddlesigns.push_back(1);
		      return 2;
		    }
		    //if(v1>=iso_val) bsaddlesigns.push_back(0); else bsaddlesigns.push_back(1);
		    return 1;
		  }
		  else 
		  {
		    if (isInCube(bsaddle2)) 
			{
		      for (int i=0;i<3;i++) temp_vtx[i]=bsaddle2[i];
		      for (int i=0;i<3;i++) bsaddle2[i]=bsaddle1[i];
		      for (int i=0;i<3;i++) bsaddle1[i]=temp_vtx[i];
			  temp_val=v1;
			  v1 = v2;
			  v2 = temp_val;
			  //if(v1>=iso_val) bsaddlesigns.push_back(0); else bsaddlesigns.push_back(1);
		      return 1;
		    }
		    return 0;
		}
	}	
}



int computeBodySaddleOld(double* bval, double* bsaddle1, double* bsaddle2)
{
  int i;
  double a,b,c,d,e,f,g,h,k0,k1,k2,zc1,zc2,z;
  double f000,f100,f010,f110,f001,f101,f011,f111;
  
  f000=bval[0];
  f100=bval[1];
  f010=bval[4];
  f110=bval[5];
  f001=bval[2];
  f101=bval[3];
  f011=bval[6];
  f111=bval[7];
  
  /*
    a= f000;
    b=-f000+f100;
    c=-f000+f010;
    d=-f000+f001;
    e=f000-f010-f100+f110;
    f=f000-f001-f100+f101;
    g=f000-f001-f010+f011;
    h=-f000+f001+f010-f011
    +f100-f101-f110+f111;
  */
  
  
  h= f000;
  e=-f000+f100;
  f=-f000+f010;
  g=-f000+f001;
  b=f000-f010-f100+f110;
  d=f000-f001-f100+f101;
  c=f000-f001-f010+f011;
  a=-f000+f001+f010-f011
    +f100-f101-f110+f111;
  
  //k0=c*f-d*e;
  //k1=g*f-d*h;
  //k2=g*e-c*h;
  
  double xc1,xc2;
  
  xc1= -c/a;
  xc2= sqrt( (a*f-b*c)*(c*d-a*g)/(a*e-b*d) ) /a;
  
  double x;
  x=bsaddle1[0]=xc1+xc2;
  bsaddle1[1]=- (b*x+f)/(a*x+c);
  bsaddle1[2]=- (d*x+g)/(a*x+c);
  
  x=bsaddle2[0]=xc1-xc2;
  bsaddle2[1]=- (b*x+f)/(a*x+c);
  bsaddle2[2]=- (d*x+g)/(a*x+c);

  /*
    
    zc1=-e/h;
    zc2=sqrt(e*e*k1*k1 - h*sqrt(k1)*(b*k2+e*k0))/h;
    z=bsaddle1[2]=zc1+zc2;
    
    bsaddle1[0]= - (c+g*z)/(e+h*z);
    bsaddle1[1]= (k0+k1*z)/k2;
    
    z=bsaddle2[2]=zc1-zc2;
    bsaddle2[0]= - (c+g*z)/(e+h*z);
    bsaddle2[1]= (k0+k1*z)/k2;
    
  */
  double temp_vtx[3];
  
  if (isInCube(bsaddle1)) {
    if (isInCube(bsaddle2)) {
      if (getTriVal(bval,bsaddle1[0],bsaddle1[1],bsaddle1[2],1) < getTriVal(bval,bsaddle2[0],bsaddle2[1],bsaddle2[2],1)) {
	for (i=0;i<3;i++) temp_vtx[i]=bsaddle2[i];
	for (i=0;i<3;i++) bsaddle2[i]=bsaddle1[i];
	for (i=0;i<3;i++) bsaddle1[i]=temp_vtx[i];
      }
      return 2;
    }
    return 1;
  } else {
    if (isInCube(bsaddle2)) {
      for (i=0;i<3;i++) temp_vtx[i]=bsaddle2[i];
      for (i=0;i<3;i++) bsaddle2[i]=bsaddle1[i];
      for (i=0;i<3;i++) bsaddle1[i]=temp_vtx[i];
      return 1;
    }
    return 0;
  }
  
}

// ret val : 0 - no ambiguity
int cubeDecomposition(double* val, Scaffold& scaffold,double isovalue)
{
  
  int i,j,k;
  int Tet[4];
  int Pyramid[5];
  double face_val[4];
  double face_saddle[6][3];
  double body_saddle1[3];
  double body_saddle2[3];
  int face_id1,face_id2;
  int ambiguous_face_idx;
  int no_bsaddle;
  double fsaddle_val;
  double temp_val,min_val,max_val;
  int fsaddle_idx;
  int big_saddle_idx,small_saddle_idx;
  
  int saddle_idx;
  int big_fsaddle,small_fsaddle;
  int big_vid,small_vid;
  int count;
  
  scaffold.clear();
  int diag_vtx;
  int face_count,next_face;
  int face_array[6],sort_face_array[6];
  int sort_first;
  int q;
  
  bool end_flag=false;
  int fi;
  int connect_no;
  bool connect_two_flag;
  bool isFaceSaddle[6];
  int no_fsaddle=0;
  int vtx_count=0;
  for (i=0;i<8;i++) scaffold.addVert(v_pos[i]);
  for (i=0;i<6;i++) {
    for (j=0;j<4;j++) {
      face_val[j]=val[face_set[i][j]];
      //for (k=0;k<3;k++) face_vpos[j][k] = v_pos[face_set[i][j]][k];
    }
    if (computeFaceSaddle(face_val,face_saddle[i],i)) {
      isFaceSaddle[i]=true;
      no_fsaddle++;
    } else isFaceSaddle[i]=false;
    scaffold.addVert(face_saddle[i]);
  }
  no_bsaddle=computeBodySaddle(val,body_saddle1,body_saddle2);
  scaffold.addVert(body_saddle1);
  scaffold.addVert(body_saddle2);
  scaffold.addVert(center_pos);
  
  double value;
  for (i=0;i<scaffold.nvert;i++) {
    if (i<8) value=val[i];
    else value=getTriVal(val,scaffold.vertices[i][0],scaffold.vertices[i][1],scaffold.vertices[i][2],1);
    scaffold.func[i]=value;
  }
  
  int ambiguous_face_no=0;
  
  int six_tet[6][4] = {
    {4,5,7,3},
    {4,5,1,3},
    {4,6,7,3},
    {4,6,2,3},
    {4,0,1,3},
    {4,0,2,3}
  };
  
  cout << no_fsaddle << " " << no_bsaddle << endl;

  //if (isType_No_Saddle(val,isovalue)) {
  if (no_fsaddle==0 && no_bsaddle==0) {
    cout << "Here arand" << endl;
    for (i=0;i<6;i++) {
      for (j=0;j<4;j++) {
	Tet[j]=six_tet[i][j];
      }
      scaffold.addTet(Tet);
    }
    //} else if (isType_SixFSaddle(val,isovalue)) {
  } else if (no_fsaddle==6) {
          
    // type 4
    int positive_fsaddle_no=0;
    big_vid=getBiggestVertex(val);
    small_vid=getSmallestVertex(val);
    
    big_fsaddle=getThirdBiggestFaceSaddle(scaffold.func,big_vid);
    small_fsaddle=getThirdSmallestFaceSaddle(scaffold.func,small_vid);
    
    
    ambiguous_face_idx=-1;
    for (i=0;i<6;i++) {
      
      for (j=0;j<4;j++) {
	face_val[j]=val[face_set[i][j]];
      }
      if (isFaceAmbiguous(face_val,isovalue)) {
	ambiguous_face_idx=i;
	ambiguous_face_no++;
	if (scaffold.func[CUBE_VTX_NO+i]>=isovalue) positive_fsaddle_no++;
      }
    }
        
    if (no_bsaddle==0) {
      
      Pyramid[0]=CUBE_VTX_NO+0;
      Pyramid[1]=CUBE_VTX_NO+1;
      Pyramid[2]=CUBE_VTX_NO+2;
      Pyramid[3]=CUBE_VTX_NO+3;
      Pyramid[4]=CUBE_VTX_NO+4;
      
      Tet[0]=Pyramid[0];
      Tet[1]=Pyramid[1];
      Tet[2]=Pyramid[2];
      Tet[3]=Pyramid[3];
      scaffold.addTet(Tet);
      
      Tet[0]=Pyramid[0];
      Tet[1]=Pyramid[1];
      Tet[2]=Pyramid[3];
      Tet[3]=Pyramid[4];
      scaffold.addTet(Tet);
      
      Pyramid[0]=CUBE_VTX_NO+5;
      Pyramid[1]=CUBE_VTX_NO+1;
      Pyramid[2]=CUBE_VTX_NO+2;
      Pyramid[3]=CUBE_VTX_NO+3;
      Pyramid[4]=CUBE_VTX_NO+4;
      
      Tet[0]=Pyramid[0];
      Tet[1]=Pyramid[1];
      Tet[2]=Pyramid[2];
      Tet[3]=Pyramid[3];
      scaffold.addTet(Tet);
      
      Tet[0]=Pyramid[0];
      Tet[1]=Pyramid[1];
      Tet[2]=Pyramid[3];
      Tet[3]=Pyramid[4];
      scaffold.addTet(Tet);
      
      count=0;
      for (i=0;i<6;i++) {
	for (j=i+1;j<6;j++) {
	  if (shared_face_edge[i][j][0]!=-1 && shared_face_edge[i][j][1]!=-1) {
	    Tet[0]= CUBE_VTX_NO+i;
	    Tet[1]= CUBE_VTX_NO+j;
	    Tet[2]= shared_face_edge[i][j][0];
	    Tet[3]= shared_face_edge[i][j][1];
	    scaffold.addTet(Tet);
	    count++;
	  }
	}
      }
      //assert(count==12);
      
      Tet[0]=CUBE_VTX_NO+0;
      Tet[1]=CUBE_VTX_NO+1;
      Tet[2]=CUBE_VTX_NO+2;
      Tet[3]=1;
      scaffold.addTet(Tet);
      
      Tet[0]=CUBE_VTX_NO+0;
      Tet[1]=CUBE_VTX_NO+2;
      Tet[2]=CUBE_VTX_NO+3;
      Tet[3]=3;
      scaffold.addTet(Tet);
      
      Tet[0]=CUBE_VTX_NO+0;
      Tet[1]=CUBE_VTX_NO+3;
      Tet[2]=CUBE_VTX_NO+4;
      Tet[3]=2;
      scaffold.addTet(Tet);
      
      Tet[0]=CUBE_VTX_NO+0;
      Tet[1]=CUBE_VTX_NO+4;
      Tet[2]=CUBE_VTX_NO+1;
      Tet[3]=0;
      scaffold.addTet(Tet);
      
      Tet[0]=CUBE_VTX_NO+5;
      Tet[1]=CUBE_VTX_NO+1;
      Tet[2]=CUBE_VTX_NO+2;
      Tet[3]=5;
      scaffold.addTet(Tet);
      
      Tet[0]=CUBE_VTX_NO+5;
      Tet[1]=CUBE_VTX_NO+2;
      Tet[2]=CUBE_VTX_NO+3;
      Tet[3]=7;
      scaffold.addTet(Tet);
      
      Tet[0]=CUBE_VTX_NO+5;
      Tet[1]=CUBE_VTX_NO+3;
      Tet[2]=CUBE_VTX_NO+4;
      Tet[3]=6;
      scaffold.addTet(Tet);
      
      Tet[0]=CUBE_VTX_NO+5;
      Tet[1]=CUBE_VTX_NO+4;
      Tet[2]=CUBE_VTX_NO+1;
      Tet[3]=4;
      scaffold.addTet(Tet);
      
      //end_flag=true;
      
      
      
    } else {
      
      getBigSmallVertex(scaffold,val,scaffold.func[BODY_SADDLE1],scaffold.func[BODY_SADDLE2],&big_vid,&small_vid);
      determine_two_face(scaffold.func,body_saddle1, body_saddle2, &face_id1 , &face_id2, big_vid, small_vid, no_bsaddle);
      
      
      if (no_bsaddle==2) {
	
	// body_saddle1 is bigger than body_saddle2
	//if (body_saddle1<body_saddle2) switch_saddle(body_saddle,1,face_id1,....);
	
	big_saddle_idx=BODY_SADDLE1;
	small_saddle_idx=BODY_SADDLE2;
      } else if (no_bsaddle==1) {
	if (scaffold.func[BODY_SADDLE2]<scaffold.func[BODY_SADDLE1]) {
	  big_saddle_idx=BODY_SADDLE1;
	  small_saddle_idx=CUBE_VTX_NO+face_id2;
	} else {
	  small_saddle_idx=BODY_SADDLE1;
	  big_saddle_idx=CUBE_VTX_NO+face_id1;
	}	
      } //else assert(0);
     
		
      if (face_id1!=big_saddle_idx-CUBE_VTX_NO) {
	for (j=0;j<4;j++) {
	  Tet[0]=big_saddle_idx;
	  Tet[1]=CUBE_VTX_NO + face_id1;
	  Tet[2]=face_set[face_id1][face_line[j][0]];
	  Tet[3]=face_set[face_id1][face_line[j][1]];
	  scaffold.addTet(Tet);
	}
      }
      
      if (face_id2!=small_saddle_idx-CUBE_VTX_NO) {
	for (j=0;j<4;j++) {
	  Tet[0]=small_saddle_idx;
	  Tet[1]=CUBE_VTX_NO + face_id2;
	  Tet[2]=face_set[face_id2][face_line[j][0]];
	  Tet[3]=face_set[face_id2][face_line[j][1]];
	  scaffold.addTet(Tet);
	}
      }
      
      
      // pick four faces
      face_count=0;
      next_face=0;
      for (i=0;i<6;i++) {
	if (i==face_id1 || i==face_id2) continue;
	face_array[face_count++]=i;
      }
      
      // sort the four faces
      sort_face_array[0]=face_array[0];
      for (i=1;i<4;i++) {
	if (is_common_edge(sort_face_array[0],face_array[i]))  {
	  next_face=face_array[i];
	  if ( is_contain(sort_face_array[0],big_vid) && is_contain(next_face,big_vid) ) {
	    sort_first=BIG;
	    sort_face_array[1]=next_face;
	    break;
	  } else if ( is_contain(sort_face_array[0],small_vid) && is_contain(next_face,small_vid) )  {
	    sort_first=SMALL;
	    sort_face_array[1]=next_face;
	    break;
	  }
	}
      }
      for (i=1;i<4;i++) {
	if (is_common_edge(sort_face_array[1],face_array[i])) {
	  next_face=face_array[i];
	  sort_face_array[2]=next_face;
	  break;
	}
      }
      
      for (i=0;i<4;i++) {
	if (face_array[i]!=sort_face_array[0] && face_array[i]!=sort_face_array[1] && face_array[i]!=sort_face_array[2]) {
	  sort_face_array[3]=face_array[i];
	  break;
	}
      }
      
      if (sort_first==BIG) {
	for (q=0;q<4;q++) {
	  i=sort_face_array[q];
	  
	  if (q==0) {
	    for (k=0;k<2;k++) {
	      Tet[0]=small_saddle_idx;
	      Tet[1]=big_vid;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(big_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	    for (k=0;k<2;k++) {
	      Tet[0]=big_saddle_idx;
	      Tet[1]=small_saddle_idx;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(big_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	    diag_vtx=face_diagonal(i,big_vid);
	    
	    for (k=0;k<2;k++) {
	      Tet[0]=big_saddle_idx;
	      Tet[1]=diag_vtx;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(big_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	  } else if (q==1) {
	    for (j=0;j<4;j++) {
	      Tet[0]=small_saddle_idx;
	      Tet[1]=CUBE_VTX_NO + i;
	      Tet[2]=face_set[i][face_line[j][0]];
	      Tet[3]=face_set[i][face_line[j][1]];
	      scaffold.addTet(Tet);
	    }
	    Tet[0]= big_saddle_idx;
	    Tet[1]= small_saddle_idx;
	    Tet[2]= shared_face_edge[i][face_id1][0];
	    Tet[3]= shared_face_edge[i][face_id1][1];
	    scaffold.addTet(Tet);
	  } else if (q==2) {
	    for (k=0;k<2;k++) {
	      Tet[0]=big_saddle_idx;
	      Tet[1]=small_vid;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(small_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	    for (k=0;k<2;k++) {
	      Tet[0]=small_saddle_idx;
	      Tet[1]=big_saddle_idx;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(small_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	    diag_vtx=face_diagonal(i,small_vid);
	    
	    for (k=0;k<2;k++) {
	      Tet[0]=small_saddle_idx;
	      Tet[1]=diag_vtx;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(small_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	  } else if (q==3) {
	    for (j=0;j<4;j++) {
	      Tet[0]=big_saddle_idx;
	      Tet[1]=CUBE_VTX_NO + i;
	      Tet[2]=face_set[i][face_line[j][0]];
	      Tet[3]=face_set[i][face_line[j][1]];
	      scaffold.addTet(Tet);
	    }
	    Tet[0]= big_saddle_idx;
	    Tet[1]= small_saddle_idx;
	    Tet[2]= shared_face_edge[i][face_id2][0];
	    Tet[3]= shared_face_edge[i][face_id2][1];
	    scaffold.addTet(Tet);
	    
	  }
	}
	
      } else {
	for (q=0;q<4;q++) {
	  i=sort_face_array[q];
	  
	  if (q==0) {
	    for (k=0;k<2;k++) {
	      Tet[0]=big_saddle_idx;
	      Tet[1]=small_vid;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(small_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	    for (k=0;k<2;k++) {
	      Tet[0]=small_saddle_idx;
	      Tet[1]=big_saddle_idx;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(small_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	    diag_vtx=face_diagonal(i,small_vid);
	    
	    for (k=0;k<2;k++) {
	      Tet[0]=small_saddle_idx;
	      Tet[1]=diag_vtx;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(small_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	  } else if (q==1) {
	    for (j=0;j<4;j++) {
	      Tet[0]=big_saddle_idx;
	      Tet[1]=CUBE_VTX_NO + i;
	      Tet[2]=face_set[i][face_line[j][0]];
	      Tet[3]=face_set[i][face_line[j][1]];
	      scaffold.addTet(Tet);
	    }
	    Tet[0]= big_saddle_idx;
	    Tet[1]= small_saddle_idx;
	    Tet[2]= shared_face_edge[i][face_id2][0];
	    Tet[3]= shared_face_edge[i][face_id2][1];
	    scaffold.addTet(Tet);
	    
	  } else if (q==2) {
	    for (k=0;k<2;k++) {
	      Tet[0]=small_saddle_idx;
	      Tet[1]=big_vid;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(big_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	    for (k=0;k<2;k++) {
	      Tet[0]=big_saddle_idx;
	      Tet[1]=small_saddle_idx;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(big_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	    diag_vtx=face_diagonal(i,big_vid);
	    
	    for (k=0;k<2;k++) {
	      Tet[0]=big_saddle_idx;
	      Tet[1]=diag_vtx;
	      Tet[2]=CUBE_VTX_NO + i;
	      Tet[3]=Vtx_Edge(big_vid,i,k);
	      scaffold.addTet(Tet);
	    }
	    
	  } else if (q==3) {
	    for (j=0;j<4;j++) {
	      Tet[0]=small_saddle_idx;
	      Tet[1]=CUBE_VTX_NO + i;
	      Tet[2]=face_set[i][face_line[j][0]];
	      Tet[3]=face_set[i][face_line[j][1]];
	      scaffold.addTet(Tet);
	    }
	    Tet[0]= big_saddle_idx;
	    Tet[1]= small_saddle_idx;
	    Tet[2]= shared_face_edge[i][face_id1][0];
	    Tet[3]= shared_face_edge[i][face_id1][1];
	    scaffold.addTet(Tet);
	  } 
	}
	
      }
    }		       
  } else if (no_bsaddle==1) {
    
    Tet[0]=BODY_SADDLE1;
    
    for (i=0;i<6;i++) {
      
      //if (isFaceAmbiguous(face_val,isovalue)) {
      if (isFaceSaddle[i]) {
	for (j=0;j<4;j++) {					
	  Tet[1]=CUBE_VTX_NO+i;
	  Tet[2]=face_set[i][face_line[j][0]];
	  Tet[3]=face_set[i][face_line[j][1]];
	  scaffold.addTet(Tet);
	}
      } else {
	for (j=0;j<2;j++) {
	  Tet[1]=sixfold_tet[i][0];
	  Tet[2]=sixfold_tet[i][1+j];
	  Tet[3]=sixfold_tet[i][2+j];
	  scaffold.addTet(Tet);
	}				
      }
      
    }
    
    //} else if (isType_NoBSaddle(val,isovalue)) {
  } else if (no_fsaddle>0) {
    
    vtx_count=0;
    double max_sval,min_sval;
    int max_fid,min_fid;
    
    if (no_fsaddle==4) {
      min_sval=10000;
      for (i=0;i<6;i++) {
	if (isFaceSaddle[i]) {
	  if (min_sval>scaffold.func[CUBE_VTX_NO+i]) {
	    min_sval=scaffold.func[CUBE_VTX_NO+i];
	    min_fid=i;
	  }
	}
      }
      max_sval=-10000;
      for (i=0;i<6;i++) {
	if (isFaceSaddle[i]) {
	  if (max_sval<scaffold.func[CUBE_VTX_NO+i]) {
	    max_sval=scaffold.func[CUBE_VTX_NO+i];
	    max_fid=i;
	  }
	}
      }
      for (i=0;i<6;i++) {
	if (isFaceSaddle[i]) {
	  if (i!=max_fid && i!=min_fid) {
	    Tet[0]=CUBE_VTX_NO+i;
	    break;
	  }
	}
      }
    } else if (no_fsaddle==3) {
      for (i=0;i<8;i++) {
	if ((val[i] >= val[vtx_connectivity[i][0]]) && 
	    (val[i] >= val[vtx_connectivity[i][1]]) && 
	    (val[i] >= val[vtx_connectivity[i][2]])) {
	  vtx_count++;
	}
      }
      
      if (vtx_count==3) { //positive=3
	
	min_sval=10000;
	for (i=0;i<6;i++) {
	  if (isFaceSaddle[i]) {
	    if (min_sval>scaffold.func[CUBE_VTX_NO+i]) {
	      min_sval=scaffold.func[CUBE_VTX_NO+i];
	      min_fid=i;
	    }
	  }
	}
	Tet[0]=CUBE_VTX_NO+min_fid;
	
      } else { // negative=3
	max_sval=-10000;
	for (i=0;i<6;i++) {
	  if (isFaceSaddle[i]) {
	    if (max_sval<scaffold.func[CUBE_VTX_NO+i]) {
	      max_sval=scaffold.func[CUBE_VTX_NO+i];
	      max_fid=i;
	    }
	  }
	}
	Tet[0]=CUBE_VTX_NO+max_fid;
	
      }
                 
    } else {
      
      for (i=0;i<6;i++) {
	if (isFaceSaddle[i]) {
	  Tet[0]=CUBE_VTX_NO+i;
	  break;
	}
      }
      
    }
    
    for (i=0;i<6;i++) {
      
      if (Tet[0]==CUBE_VTX_NO+i) continue;
      //if (isFaceAmbiguous(face_val,isovalue)) {
      if (isFaceSaddle[i]) {
	for (j=0;j<4;j++) {					
	  Tet[1]=CUBE_VTX_NO+i;
	  Tet[2]=face_set[i][face_line[j][0]];
	  Tet[3]=face_set[i][face_line[j][1]];
	  scaffold.addTet(Tet);
	}
      } else {
	for (j=0;j<2;j++) {
	  Tet[1]=sixfold_tet[i][0];
	  Tet[2]=sixfold_tet[i][1+j];
	  Tet[3]=sixfold_tet[i][2+j];
	  scaffold.addTet(Tet);
	}				
      }
      
    }
    
    
  } //else assert(0); 
  
  // big comment starts here...
  /*
  if (isType1(val,isovalue)) {
    return 0;
  } else if (isType2(val,isovalue,no_bsaddle)) {    
    ambiguous_face_idx=-1;
    for (i=0;i<6;i++) {
      
      for (j=0;j<4;j++) {
	face_val[j]=val[face_set[i][j]];
      }
      if (isFaceAmbiguous(face_val,isovalue)) {
	ambiguous_face_idx=i;
	ambiguous_face_no++;
      }
    }
    
    if (ambiguous_face_no==0) {
      //assert(no_bsaddle==1);
      if (getPositiveNo(val,isovalue)==2) {
	if (vertex_sign(scaffold.func[BODY_SADDLE1],isovalue)==SIGN_POSITIVE) {
	  Tet[0]=BODY_SADDLE1;
	  Pyramid[0]=BODY_SADDLE1;
	} else {
	  return 0;
	}
      } else if (getPositiveNo(val,isovalue)==6) {
	if (vertex_sign(scaffold.func[BODY_SADDLE1],isovalue)==SIGN_NEGATIVE) {
	  Tet[0]=BODY_SADDLE1;
	  Pyramid[0]=BODY_SADDLE1;
	} else {
	  return 0;
	}
	
      } //else assert(0);
    } else if (ambiguous_face_no==1) {
      if (getPositiveNo(val,isovalue)==2 || getPositiveNo(val,isovalue)==6) {
	// case 3
	Tet[0]=CUBE_VTX_NO+ambiguous_face_idx;
	Pyramid[0]=CUBE_VTX_NO+ambiguous_face_idx;
      } else if (getPositiveNo(val,isovalue)==3) {
	// case 6_positive
	if (no_bsaddle==0) {
	  Tet[0]=CUBE_VTX_NO+ambiguous_face_idx;
	  Pyramid[0]=CUBE_VTX_NO+ambiguous_face_idx;
	} else {
	  if (vertex_sign(scaffold.func[BODY_SADDLE1],isovalue)==SIGN_POSITIVE) {
	    if (vertex_sign(scaffold.func[CUBE_VTX_NO+ambiguous_face_idx],isovalue)==SIGN_POSITIVE) {
	      Tet[0]=CUBE_VTX_NO+ambiguous_face_idx;
	      Pyramid[0]=CUBE_VTX_NO+ambiguous_face_idx;
	    } else {
	      Tet[0]=BODY_SADDLE1;
	      Pyramid[0]=BODY_SADDLE1;
	    }
	  } else {
	    Tet[0]=CUBE_VTX_NO+ambiguous_face_idx;
	    Pyramid[0]=CUBE_VTX_NO+ambiguous_face_idx;
	  }
	}
      } else if (getPositiveNo(val,isovalue)==5) {
	// case 6_negative
	if (no_bsaddle==0) {
	  Tet[0]=CUBE_VTX_NO+ambiguous_face_idx;
	  Pyramid[0]=CUBE_VTX_NO+ambiguous_face_idx;
	} else {
	  if (vertex_sign(scaffold.func[BODY_SADDLE1],isovalue)==SIGN_NEGATIVE) {
	    if (vertex_sign(scaffold.func[CUBE_VTX_NO+ambiguous_face_idx],isovalue)==SIGN_NEGATIVE) {
	      Tet[0]=CUBE_VTX_NO+ambiguous_face_idx;
	      Pyramid[0]=CUBE_VTX_NO+ambiguous_face_idx;
	    } else {
	      Tet[0]=BODY_SADDLE1;
	      Pyramid[0]=BODY_SADDLE1;
	    }
	  } else {
	    Tet[0]=CUBE_VTX_NO+ambiguous_face_idx;
	    Pyramid[0]=CUBE_VTX_NO+ambiguous_face_idx;
	  }
	}
      }
      
    } else if (ambiguous_face_no==2) {
      if (no_bsaddle==0) {
	Tet[0]=CUBE_VTX_NO+ambiguous_face_idx;
	Pyramid[0]=CUBE_VTX_NO+ambiguous_face_idx;
      } else if (no_bsaddle==1) {
	for (i=0;i<6;i++) {
	  for (j=0;j<4;j++) {
	    face_val[j]=val[face_set[i][j]];
	  }
	  
	  if (isFaceAmbiguous(face_val,isovalue)) {
	    if (vertex_sign(scaffold.func[BODY_SADDLE1],isovalue) == vertex_sign(scaffold.func[CUBE_VTX_NO+i],isovalue)) {
	      break;
	    }
	  }
	}
	if (i==6) {
	  Tet[0]=BODY_SADDLE1;
	  Pyramid[0]=BODY_SADDLE1;
	} else {
	  Tet[0]=CUBE_VTX_NO+ambiguous_face_idx;
	  Pyramid[0]=CUBE_VTX_NO+ambiguous_face_idx;
	}
	
      } //else assert(0);
    }   
    
    for (i=0;i<6;i++) {
      for (j=0;j<4;j++) {
	face_val[j]=val[face_set[i][j]];
      }
      if (isFaceAmbiguous(face_val,isovalue)) {
	
	for (j=0;j<4;j++) {					
	  Tet[1]=CUBE_VTX_NO+i;
	  Tet[2]=face_set[i][face_line[j][0]];
	  Tet[3]=face_set[i][face_line[j][1]];
	  scaffold.addTet(Tet);
	}	
	
      } else {
	for (j=0;j<4;j++) Pyramid[j+1]=face_set[i][j];
				scaffold.addPyramid(Pyramid);
      }
    }
      
  } else if (isType3(val,isovalue)) { // type 3

    min_val=1000000000;
    max_val=-1000000000;
    switch (no_bsaddle) {
    case 0 :
      if (getPositiveNo(val,isovalue)==3) {
	for (i=0;i<6;i++) {
	  for (j=0;j<4;j++) {
	    face_val[j]=val[face_set[i][j]];
	  }
	  if (isFaceAmbiguous(face_val,isovalue)) {
	    temp_val=scaffold.func[CUBE_VTX_NO+i];
	    if (temp_val<min_val) {
	      min_val=temp_val;
	      saddle_idx=i;
	    }
	  }
	}
	
      } else if (getPositiveNo(val,isovalue)==5) {
	for (i=0;i<6;i++) {
	  for (j=0;j<4;j++) {
	    face_val[j]=val[face_set[i][j]];
	  }
	  if (isFaceAmbiguous(face_val,isovalue)) {
	    temp_val=scaffold.func[CUBE_VTX_NO+i];
	    if (temp_val>max_val) {
	      max_val=temp_val;
	      saddle_idx=i;
	    }
	  }
	}       	
	
      } //else assert(0);
                 
      for (i=0;i<6;i++) {
	
	for (j=0;j<4;j++) {
	  face_val[j]=val[face_set[i][j]];
	}
	
	if (isFaceAmbiguous(face_val,isovalue)) {
	  
	  if (saddle_idx==i) {
	    continue;						
	  } else {
	    for (j=0;j<4;j++) {					
	      Tet[0]=CUBE_VTX_NO+saddle_idx;
	      Tet[1]=CUBE_VTX_NO+i;
	      Tet[2]=face_set[i][face_line[j][0]];
	      Tet[3]=face_set[i][face_line[j][1]];
	      scaffold.addTet(Tet);
	    }
	  }				
	  
	} else {
	  Pyramid[0]=CUBE_VTX_NO+saddle_idx;
	  for (j=0;j<4;j++) Pyramid[j+1]=face_set[i][j];
	  scaffold.addPyramid(Pyramid);
	}
      }
      break;
    case 1 :
      if (getPositiveNo(val,isovalue)==3) {
	for (i=0;i<6;i++) {
	  for (j=0;j<4;j++) {
	    face_val[j]=val[face_set[i][j]];
	  }
	  if (isFaceAmbiguous(face_val,isovalue)) {
	    temp_val=scaffold.func[CUBE_VTX_NO+i];
	    if (temp_val<min_val) {
	      min_val=temp_val;
	      fsaddle_idx=i;
	    }
	  }
	}
	
      } else if (getPositiveNo(val,isovalue)==5) {
	for (i=0;i<6;i++) {
	  for (j=0;j<4;j++) {
	    face_val[j]=val[face_set[i][j]];
	  }
	  if (isFaceAmbiguous(face_val,isovalue)) {
	    temp_val=scaffold.func[CUBE_VTX_NO+i];
	    if (temp_val>max_val) {
	      max_val=temp_val;
	      fsaddle_idx=i;
	    }
	  }
	}
      } //else assert(0);
      
      if (getPositiveNo(val,isovalue) == 3) {
	if (scaffold.func[BODY_SADDLE1]>scaffold.func[BODY_SADDLE2]) saddle_idx=BODY_SADDLE1;
	else saddle_idx=CUBE_VTX_NO+fsaddle_idx;
	
      } else {
	if (scaffold.func[BODY_SADDLE1]<scaffold.func[BODY_SADDLE2]) saddle_idx=BODY_SADDLE1;
	else saddle_idx=CUBE_VTX_NO+fsaddle_idx;
      }
      
      if (saddle_idx==BODY_SADDLE1) {
	for (i=0;i<6;i++) {
	  for (j=0;j<4;j++) {
	    face_val[j]=val[face_set[i][j]];
	  }
	  if (isFaceAmbiguous(face_val,isovalue)) {
	    if (vertex_sign(scaffold.func[BODY_SADDLE1],isovalue)==vertex_sign(scaffold.func[CUBE_VTX_NO+i],isovalue))
	      break;
	  }
	}
	if (i!=6) saddle_idx=CUBE_VTX_NO+fsaddle_idx;
	
      }
      
      for (i=0;i<6;i++) {
	for (j=0;j<4;j++) {
	  face_val[j]=val[face_set[i][j]];
	}
	
	if (saddle_idx==i+CUBE_VTX_NO) continue;
	
	if (isFaceAmbiguous(face_val,isovalue)) {
	  for (j=0;j<4;j++) {					
	    Tet[0]=saddle_idx;
	    Tet[1]=CUBE_VTX_NO+i;
	    Tet[2]=face_set[i][face_line[j][0]];
	    Tet[3]=face_set[i][face_line[j][1]];
	    scaffold.addTet(Tet);
	  }				
	  
	} else {
	  Pyramid[0]=saddle_idx;
	  for (j=0;j<4;j++) Pyramid[j+1]=face_set[i][j];
	  scaffold.addPyramid(Pyramid);
	}
      }
      break;
    case 2 :
      if (getPositiveNo(val,isovalue)==3) {
	if (scaffold.func[BODY_SADDLE1]>scaffold.func[BODY_SADDLE2]) {
	  saddle_idx=BODY_SADDLE1;
	} else saddle_idx=BODY_SADDLE2;
      } else if (getPositiveNo(val,isovalue)==5) {
	if (scaffold.func[BODY_SADDLE1]>scaffold.func[BODY_SADDLE2]) {
	  saddle_idx=BODY_SADDLE2;
	} else saddle_idx=BODY_SADDLE1;
      } //else assert(0);
      
      
      if (saddle_idx==BODY_SADDLE1 || saddle_idx==BODY_SADDLE2) {
	for (i=0;i<6;i++) {
	  for (j=0;j<4;j++) {
	    face_val[j]=val[face_set[i][j]];
	  }
	  if (isFaceAmbiguous(face_val,isovalue)) {
	    if (vertex_sign(scaffold.func[saddle_idx],isovalue)==vertex_sign(scaffold.func[CUBE_VTX_NO+i],isovalue))
	      break;
	  }
	}
	if (i!=6) saddle_idx=CUBE_VTX_NO+fsaddle_idx;
	
      }
      
      
      
      for (i=0;i<6;i++) {
	for (j=0;j<4;j++) {
	  face_val[j]=val[face_set[i][j]];
	}
	
	if (isFaceAmbiguous(face_val,isovalue)) {
	  for (j=0;j<4;j++) {					
	    Tet[0]=saddle_idx;
	    Tet[1]=CUBE_VTX_NO+i;
	    Tet[2]=face_set[i][face_line[j][0]];
	    Tet[3]=face_set[i][face_line[j][1]];
	    scaffold.addTet(Tet);
	  }				
	  
	} else {
	  Pyramid[0]=saddle_idx;
	  for (j=0;j<4;j++) Pyramid[j+1]=face_set[i][j];
	  scaffold.addPyramid(Pyramid);
	}
      }
      break;
    default : 1;//assert(0);
    }
       
  } else { 
    // type 4
    int positive_fsaddle_no=0;
    big_vid=getBiggestVertex(val);
    small_vid=getSmallestVertex(val);
    
    big_fsaddle=getThirdBiggestFaceSaddle(scaffold.func,big_vid);
    small_fsaddle=getThirdSmallestFaceSaddle(scaffold.func,small_vid);
        
    ambiguous_face_idx=-1;
    for (i=0;i<6;i++) {
      
      for (j=0;j<4;j++) {
	face_val[j]=val[face_set[i][j]];
      }
      if (isFaceAmbiguous(face_val,isovalue)) {
	ambiguous_face_idx=i;
	ambiguous_face_no++;
	if (scaffold.func[CUBE_VTX_NO+i]>=isovalue) positive_fsaddle_no++;
      }
    }
    
    if (positive_fsaddle_no==3) {
      
      if (no_bsaddle==2) {

	// arand commented... this was broken...
	//determine_two_face(body_saddle1, body_saddle2, &face_id1 , &face_id2);

	// body_saddle1 is bigger than body_saddle2
	//if (body_saddle1<body_saddle2) switch_saddle(body_saddle,1,face_id1,....);
	
	big_saddle_idx=BODY_SADDLE1;
	small_saddle_idx=BODY_SADDLE2;
	
	getBigSmallVertex(scaffold,val,scaffold.func[BODY_SADDLE1],scaffold.func[BODY_SADDLE2],&big_vid,&small_vid);
	
	
	
	if (face_id1!=big_saddle_idx-CUBE_VTX_NO) {
	  for (j=0;j<4;j++) {
	    Tet[0]=big_saddle_idx;
	    Tet[1]=CUBE_VTX_NO + face_id1;
	    Tet[2]=face_set[face_id1][face_line[j][0]];
	    Tet[3]=face_set[face_id1][face_line[j][1]];
	    scaffold.addTet(Tet);
	  }
	}
	
	if (face_id2!=small_saddle_idx-CUBE_VTX_NO) {
	  for (j=0;j<4;j++) {
	    Tet[0]=small_saddle_idx;
	    Tet[1]=CUBE_VTX_NO + face_id2;
	    Tet[2]=face_set[face_id2][face_line[j][0]];
	    Tet[3]=face_set[face_id2][face_line[j][1]];
	    scaffold.addTet(Tet);
	  }
	}
	
	
	// pick four faces
	face_count=0;
	next_face=0;
	for (i=0;i<6;i++) {
	  if (i==face_id1 || i==face_id2) continue;
	  face_array[face_count++]=i;
	}
	
	// sort the four faces
	sort_face_array[0]=face_array[0];
	for (i=1;i<4;i++) {
	  if (is_common_edge(sort_face_array[0],face_array[i]))  {
	    next_face=face_array[i];
	    if ( is_contain(sort_face_array[0],big_vid) && is_contain(next_face,big_vid) ) {
	      sort_first=BIG;
	      sort_face_array[1]=next_face;
	      break;
	    } else if ( is_contain(sort_face_array[0],small_vid) && is_contain(next_face,small_vid) )  {
	      sort_first=SMALL;
	      sort_face_array[1]=next_face;
	      break;
	    }
	  }
	}
	for (i=1;i<4;i++) {
	  if (is_common_edge(sort_face_array[1],face_array[i])) {
	    next_face=face_array[i];
	    sort_face_array[2]=next_face;
	    break;
	  }
	}
	
	for (i=0;i<4;i++) {
	  if (face_array[i]!=sort_face_array[0] && face_array[i]!=sort_face_array[1] && face_array[i]!=sort_face_array[2]) {
	    sort_face_array[3]=face_array[i];
	    break;
	  }
	}
	
	if (sort_first==BIG) {
	  for (q=0;q<4;q++) {
	    i=sort_face_array[q];
	    
	    if (q==0) {
	      for (k=0;k<2;k++) {
		Tet[0]=small_saddle_idx;
		Tet[1]=big_vid;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(big_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	      for (k=0;k<2;k++) {
		Tet[0]=big_saddle_idx;
		Tet[1]=small_saddle_idx;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(big_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	      diag_vtx=face_diagonal(i,big_vid);
	      
	      for (k=0;k<2;k++) {
		Tet[0]=big_saddle_idx;
		Tet[1]=diag_vtx;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(big_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	    } else if (q==1) {
	      for (j=0;j<4;j++) {
		Tet[0]=small_saddle_idx;
		Tet[1]=CUBE_VTX_NO + i;
		Tet[2]=face_set[i][face_line[j][0]];
		Tet[3]=face_set[i][face_line[j][1]];
		scaffold.addTet(Tet);
	      }
	      Tet[0]= big_saddle_idx;
	      Tet[1]= small_saddle_idx;
	      Tet[2]= shared_face_edge[i][face_id1][0];
	      Tet[3]= shared_face_edge[i][face_id1][1];
	      scaffold.addTet(Tet);
	    } else if (q==2) {
	      for (k=0;k<2;k++) {
		Tet[0]=big_saddle_idx;
		Tet[1]=small_vid;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(small_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	      for (k=0;k<2;k++) {
		Tet[0]=small_saddle_idx;
		Tet[1]=big_saddle_idx;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(small_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	      diag_vtx=face_diagonal(i,small_vid);
	      
	      for (k=0;k<2;k++) {
		Tet[0]=small_saddle_idx;
		Tet[1]=diag_vtx;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(small_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	    } else if (q==3) {
	      for (j=0;j<4;j++) {
		Tet[0]=big_saddle_idx;
		Tet[1]=CUBE_VTX_NO + i;
		Tet[2]=face_set[i][face_line[j][0]];
		Tet[3]=face_set[i][face_line[j][1]];
		scaffold.addTet(Tet);
	      }
	      Tet[0]= big_saddle_idx;
	      Tet[1]= small_saddle_idx;
	      Tet[2]= shared_face_edge[i][face_id2][0];
	      Tet[3]= shared_face_edge[i][face_id2][1];
	      scaffold.addTet(Tet);
	      
	    }
	  }
	  
	} else {
	  for (q=0;q<4;q++) {
	    i=sort_face_array[q];
	    
	    if (q==0) {
	      for (k=0;k<2;k++) {
		Tet[0]=big_saddle_idx;
		Tet[1]=small_vid;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(small_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	      for (k=0;k<2;k++) {
		Tet[0]=small_saddle_idx;
		Tet[1]=big_saddle_idx;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(small_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	      diag_vtx=face_diagonal(i,small_vid);
	      
	      for (k=0;k<2;k++) {
		Tet[0]=small_saddle_idx;
		Tet[1]=diag_vtx;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(small_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	    } else if (q==1) {
	      for (j=0;j<4;j++) {
		Tet[0]=big_saddle_idx;
		Tet[1]=CUBE_VTX_NO + i;
		Tet[2]=face_set[i][face_line[j][0]];
		Tet[3]=face_set[i][face_line[j][1]];
		scaffold.addTet(Tet);
	      }
	      Tet[0]= big_saddle_idx;
	      Tet[1]= small_saddle_idx;
	      Tet[2]= shared_face_edge[i][face_id2][0];
	      Tet[3]= shared_face_edge[i][face_id2][1];
	      scaffold.addTet(Tet);
	      
	    } else if (q==2) {
	      for (k=0;k<2;k++) {
		Tet[0]=small_saddle_idx;
		Tet[1]=big_vid;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(big_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	      for (k=0;k<2;k++) {
		Tet[0]=big_saddle_idx;
		Tet[1]=small_saddle_idx;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(big_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	      diag_vtx=face_diagonal(i,big_vid);
	      
	      for (k=0;k<2;k++) {
		Tet[0]=big_saddle_idx;
		Tet[1]=diag_vtx;
		Tet[2]=CUBE_VTX_NO + i;
		Tet[3]=Vtx_Edge(big_vid,i,k);
		scaffold.addTet(Tet);
	      }
	      
	    } else if (q==3) {
	      for (j=0;j<4;j++) {
		Tet[0]=small_saddle_idx;
		Tet[1]=CUBE_VTX_NO + i;
		Tet[2]=face_set[i][face_line[j][0]];
		Tet[3]=face_set[i][face_line[j][1]];
		scaffold.addTet(Tet);
	      }
	      Tet[0]= big_saddle_idx;
	      Tet[1]= small_saddle_idx;
	      Tet[2]= shared_face_edge[i][face_id1][0];
	      Tet[3]= shared_face_edge[i][face_id1][1];
	      scaffold.addTet(Tet);
	    } 
	  }
	  
	}
	
      }
    }    
  }		

  */


  for (i=0;i<scaffold.nvert;i++) {
    for (j=0;j<3;j++) scaffold.vertices[i][j]-=0.5;
  }
  return 1; 
}


void cube_triangulate(double* val, double isovalue)//, geoframe& geom)
{

  // arand: regular marching cubes here... fix
  /*
  Contouring* onecell;
  int i;
  
  
  int Dim[3];
  double Orig[3],Span[3];
  
  for (i=0;i<3;i++) {
    Dim[i]=2;
    Span[i]=1.0f;
    Orig[i]=-.5f;
  }
  onecell=new Contouring(val,Dim,Orig,Span);
  
  onecell->getIsosurfaceGeom(isovalue,geom);
    
  delete onecell;
  */
}


void tet_edge2pos(int eid, double* vpos0,double  vertices[][3],double* func, double isovalue)
{
  int i;
  int v1,v2;
  double f1,f2,ratio;
  
  v1=tet_edge[eid][0];
  v2=tet_edge[eid][1];
  f1=func[v1];
  f2=func[v2];
  
  if (f1!=f2)
    ratio=(isovalue-f1)/(f2-f1);
  else ratio=0.5;
  
  for (i=0;i<3;i++) vpos0[i]=vertices[v1][i]*(1-ratio)+ratio*vertices[v2][i];
}

void pyramid_edge2pos(int eid, double* vpos0,double  vertices[][3],double* func, double isovalue)
{
  int i;
  int v1,v2;
  double f1,f2;
  double ratio;
  
  v1=pyramid_edge[eid][0];
  v2=pyramid_edge[eid][1];
  f1=func[v1];
  f2=func[v2];
  
  if (f1!=f2)
    ratio=(isovalue-f1)/(f2-f1);
  else ratio=0.5;
  
  for (i=0;i<3;i++) vpos0[i]=vertices[v1][i]*(1-ratio)+ratio*vertices[v2][i];
}

void tet_edge2pos_tri(int eid, double* vpos0,double  vertices[][3],double* cube_val, double isovalue)
{
  int i,v1,v2,iterate_i;
  double f1,f2;
  double vpos1[3],vpos2[3],center_pos[3];
  double center_val;
  
  v1=tet_edge[eid][0];
  v2=tet_edge[eid][1];
  
  for (i=0;i<3;i++) {
    vpos1[i]=vertices[v1][i];
    vpos2[i]=vertices[v2][i];
  }
  
  for (iterate_i=0;iterate_i<10;iterate_i++) {
    f1 = getTriVal(cube_val,vpos1[0]+0.5,vpos1[1]+0.5,vpos1[2]+0.5,1);
    f2 = getTriVal(cube_val,vpos2[0]+0.5,vpos2[1]+0.5,vpos2[2]+0.5,1);
    
    for (i=0;i<3;i++) center_pos[i] = 0.5*vpos1[i]+0.5*vpos2[i];
    
    center_val = getTriVal(cube_val, center_pos[0]+0.5  , center_pos[1]+0.5  , center_pos[2]+0.5  , 1);
    
    if (f1<isovalue) {
      if (center_val<isovalue) {
	for (i=0;i<3;i++) vpos1[i]=center_pos[i];
      } else {
	for (i=0;i<3;i++) vpos2[i] = center_pos[i];
      }
    } else if (f1>=isovalue) {
      if (center_val<isovalue) {
	for (i=0;i<3;i++) vpos2[i] = center_pos[i];
      } else {
	for (i=0;i<3;i++) vpos1[i] = center_pos[i];
      }
    }
  }
  for (i=0;i<3;i++) vpos0[i]=center_pos[i];
}


void pyramid_edge2pos_tri(int eid, double* vpos0,double  vertices[][3],double* cube_val, double isovalue)
{
  int i,v1,v2,iterate_i;
  double f1,f2;
  double vpos1[3],vpos2[3],center_pos[3];
  double center_val;

  v1=pyramid_edge[eid][0];
  v2=pyramid_edge[eid][1];
  
  for (i=0;i<3;i++) {
    vpos1[i]=vertices[v1][i];
    vpos2[i]=vertices[v2][i];
  }
  
  for (iterate_i=0;iterate_i<10;iterate_i++) {
    f1 = getTriVal(cube_val,vpos1[0]+0.5,vpos1[1]+0.5,vpos1[2]+0.5,1);
    f2 = getTriVal(cube_val,vpos2[0]+0.5,vpos2[1]+0.5,vpos2[2]+0.5,1);
    
    for (i=0;i<3;i++) center_pos[i] = 0.5*vpos1[i]+0.5*vpos2[i];
    
    center_val = getTriVal(cube_val, center_pos[0]+0.5  , center_pos[1]+0.5  , center_pos[2]+0.5  , 1);
    
    if (f1<isovalue) {
      if (center_val<isovalue) {
	for (i=0;i<3;i++) vpos1[i]=center_pos[i];
      } else {
	for (i=0;i<3;i++) vpos2[i] = center_pos[i];
      }
    } else if (f1>=isovalue) {
      if (center_val<isovalue) {
	for (i=0;i<3;i++) vpos2[i] = center_pos[i];
      } else {
	for (i=0;i<3;i++) vpos1[i] = center_pos[i];
      }
    }
  }
  for (i=0;i<3;i++) vpos0[i]=center_pos[i];
}



int smooth_flag=1;
void tet_triangulate(double  vertices[][3],Cell* c,double * corner, int tet[][4],int ntet , double isovalue) //, geoframe& gf)
{

  cout << "NTET " << ntet << endl;
  
  double * func = c->values;

  double val[4];
  int code;
  int e0,e1,e2;
  double vpos0[3],vpos1[3],vpos2[3];
  int v0,v1,v2;
  int i,j,k;
  int tri_no;
  double norm[3];
  double tet_vertices[4][3];
  
  for (i=0;i<ntet;i++) {
    //for (i=0;i<2;i++) {
    for (j=0;j<4;j++) {
      val[j]=func[tet[i][j]];
      for (k=0;k<3;k++) tet_vertices[j][k]=vertices[tet[i][j]][k];
    }
    
    
    code=0;
    if (val[0]>isovalue) code|= 0x01;
    if (val[1]>isovalue) code|= 0x02;
    if (val[2]>isovalue) code|= 0x04;
    if (val[3]>isovalue) code|= 0x08;
    
    tri_no=tet_triangle[code][0];
    
    for (j=0;j<tri_no;j++) {
      
      e0=tet_triangle[code][1+j*3+0];
      e1=tet_triangle[code][1+j*3+1];
      e2=tet_triangle[code][1+j*3+2];
      
      if (smooth_flag==0) {
	tet_edge2pos(e0,vpos0,tet_vertices,val,isovalue);
	tet_edge2pos(e1,vpos1,tet_vertices,val,isovalue);
	tet_edge2pos(e2,vpos2,tet_vertices,val,isovalue);
      } else {
	tet_edge2pos_tri(e0,vpos0,tet_vertices,func,isovalue);
	tet_edge2pos_tri(e1,vpos1,tet_vertices,func,isovalue);
	tet_edge2pos_tri(e2,vpos2,tet_vertices,func,isovalue);
      }
           

      tri0.push_back(mvx.size());

      mvx.push_back(corner[0] + c->span*(1+vpos0[0]));
      mvy.push_back(corner[1] + c->span*(1+vpos0[1]));
      mvz.push_back(corner[2] + c->span*(1+vpos0[2]));

      tri1.push_back(mvx.size());
      mvx.push_back(corner[0] + c->span*(1+vpos1[0]));
      mvy.push_back(corner[1] + c->span*(1+vpos1[1]));
      mvz.push_back(corner[2] + c->span*(1+vpos1[2]));

      tri2.push_back(mvx.size());
      mvx.push_back(corner[0] + c->span*(1+vpos2[0]));
      mvy.push_back(corner[1] + c->span*(1+vpos2[1]));
      mvz.push_back(corner[2] + c->span*(1+vpos2[2]));

      /*
      cout << "(tet) NewTri: (" << vpos0[0] << ","  
	   << vpos0[1] << ","  << vpos0[2] << ")"   << endl;

      cout << "(tet) NewTri: (" << vpos1[0] << ","  
	   << vpos1[1] << ","  << vpos1[2] << ")"   << endl;

      cout << "(tet) NewTri: (" << vpos2[0] << ","  
	   << vpos2[1] << ","  << vpos2[2] << ")"   << endl;
      */

      /*
      v0=gf.AddVert(vpos0,norm);
      v1=gf.AddVert(vpos1,norm);
      v2=gf.AddVert(vpos2,norm);

      gf.AddTri(v0,v1,v2);
      */
    }
    
  }
  
}


void pyramid_triangulate(double vertices[][3], double* func, int  pyramid[][5],int npyramid, double isovalue) //,geoframe& gf)
{
  double val[5];
  int code;
  int e0,e1,e2;
  double vpos0[3],vpos1[3],vpos2[3];
  int v0,v1,v2;
  int i,j,k;
  int pyramid_no;
  double pyramid_vertices[5][3];
  double norm[3];
    
  for (i=0;i<npyramid;i++) {
    //for (i=0;i<1;i++) {
    for (j=0;j<5;j++) {
      val[j]=func[pyramid[i][j]];
      for (k=0;k<3;k++) pyramid_vertices[j][k]=vertices[pyramid[i][j]][k];
    }
    
    code=0;
    if (val[0]>isovalue) code|= 0x10;
    if (val[1]>isovalue) code|= 0x01;
    if (val[2]>isovalue) code|= 0x02;
    if (val[3]>isovalue) code|= 0x04;
    if (val[4]>isovalue) code|= 0x08;
    
    pyramid_no=pyramid_triangle[code][0];
    
    for (j=0;j<pyramid_no;j++) {
      
      e0=pyramid_triangle[code][1+j*3+0];
      e1=pyramid_triangle[code][1+j*3+1];
      e2=pyramid_triangle[code][1+j*3+2];
      
      if (smooth_flag==0) {
	pyramid_edge2pos(e0,vpos0,pyramid_vertices,val,isovalue);
	pyramid_edge2pos(e1,vpos1,pyramid_vertices,val,isovalue);
	pyramid_edge2pos(e2,vpos2,pyramid_vertices,val,isovalue);
      } else {
	pyramid_edge2pos_tri(e0,vpos0,pyramid_vertices,func,isovalue);
	pyramid_edge2pos_tri(e1,vpos1,pyramid_vertices,func,isovalue);
	pyramid_edge2pos_tri(e2,vpos2,pyramid_vertices,func,isovalue);
      }

      cout << "(pyr) NewTri: (" << vpos0[0] << ","  
	   << vpos0[1] << ","  << vpos0[2] << ")"   << endl;
      
      /* // fix this
	 v0=gf.AddVert(vpos0,norm);
      v1=gf.AddVert(vpos1,norm);
      v2=gf.AddVert(vpos2,norm);
      
      gf.AddTri(v0,v1,v2);
      */
    }
    
  }
}




int ImprovedMC(Cell* cell, double * corner, double isoval) {

  Scaffold scaf;   

  double * val = cell->values;

  int type_id=cubeDecomposition(val, scaf,isoval);
  cout << "TI " << type_id << endl;
  if (type_id==0) {
    cube_triangulate(val,isoval);//,g_frames[0]);
  } else {
    tet_triangulate(scaf.vertices,cell,corner,scaf.tet,scaf.ntet , isoval);//, g_frames[0]);
    pyramid_triangulate(scaf.vertices,scaf.func,scaf.pyramid,scaf.npyramid , isoval);//, g_frames[0]);
  }

  cout << val[0] << " " << val[1] << " " << val[2] << " " << val[3] << " " << val[4] << " " << val[5] << " " << val[6] << " " << val[7] << " " << endl;  
  cout << "Size: " << mvx.size() << endl;

  return type_id;

}
