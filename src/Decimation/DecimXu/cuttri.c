#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <Decimation/geometry.h>

int CutTriangleByTwoPlane(p1,p2,p3,n1,n2,n3,a,b,tri1)
float p1[],p2[],p3[],n1[],n2[],n3[],a,b;
Triangle tri1;
{
int   i;

if ( (p1[0] <= a) && (p2[0] <= a) && (p3[0] <= a) &&
         (p1[2] >= b) && (p2[2] >= b) && (p3[2] >= b) ) {
   return(0);
} else {
for (i = 0; i < 3; i++) {
      tri1->Vertex1[i] = p1[i];
      tri1->Vertex2[i] = p2[i];
      tri1->Vertex3[i] = p3[i];
      tri1->Normal1[i] = n1[i];
      tri1->Normal2[i] = n2[i];
      tri1->Normal3[i] = n3[i];
   }
   return(1);
}

}

void TransferFile_raw()
{
int           numbtris, i,j,k,l,index[3],numbpts,*triangles,eqal;
double        x,y,z,nx,ny,nz;
float         *points,dis; 
FILE          *fppts1, *fppts2; 


fppts1 = fopen("../../data/Triangles/tmp2","r");

fscanf(fppts1,"%d\n", &numbtris);

points = (float *) malloc(18*numbtris*sizeof(float));
triangles = (int *) malloc(3*numbtris*sizeof(int));

numbpts = 0;
for (i = 0 ; i < numbtris; i++) {
    for (l = 0; l < 3; l++) { 
       fscanf(fppts1,"%lf %lf %lf %lf %lf %lf\n",
                                &x, &y, &z, &nx, &ny, &nz);
       eqal = 0;
       for (j = 0; j < numbpts; j++) {
           k = 6*j;
           dis = fabs(points[k]     - x) + 
                 fabs(points[k + 1] - y) +
                 fabs(points[k + 2] - z); 
           if (dis < 0.00001) {
              eqal = 1;
              index[l] = j;
              break;
           }
       }
       if (eqal == 0) {
          points[6*numbpts] = x;
          points[6*numbpts + 1] = y;
          points[6*numbpts + 2] = z;
          points[6*numbpts + 3] = nx;
          points[6*numbpts + 4] = ny;
          points[6*numbpts + 5] = nz;
          index[l] = numbpts;
          numbpts = numbpts + 1;
       }
    }
    triangles[3*i]   = index[0];
    triangles[3*i+1] = index[1];
    triangles[3*i+2] = index[2];
}
fclose(fppts1);

fppts2 = fopen("../../data/Triangles/tmp3.glx","w");
fprintf(fppts2,"%10d %10d\n",numbpts, numbtris);
/*
printf("save with normal\n"); 
for (i = 0; i < numbpts; i++) {
    fprintf(fppts2, "%lf %lf %lf %lf %lf %lf\n", points[6*i], points[6*i+1], points[6*i+2],
                   points[6*i+3], points[6*i+4], points[6*i+5]);
}
*/
printf("save without normal\n"); 
for (i = 0; i < numbpts; i++) {
    fprintf(fppts2, "%lf %lf %lf\n", points[6*i], points[6*i+1], points[6*i+2]);
}

for (i = 0; i < numbtris; i++) {
    fprintf(fppts2, "%d %d %d\n", triangles[3*i],triangles[3*i+1],triangles[3*i+2]);
}

free(points); 
free(triangles);
fclose(fppts2);
}

void TransferToFat()
{
FILE          *fppts1, *fppts2;
int           i,numbpts, numbtris,ii,jj,kk;
float         det;
double        x,y,z,nx,ny,nz;


fppts1 = fopen("../../data/Triangles/tmp3.glx","r");
/*
fppts1 = fopen("../../data/Triangles/tmp","r");
*/
fppts2 = fopen("../../data/FatTri/tmp.fat","w");

fscanf(fppts1,"%d %d\n", &numbpts, &numbtris);
fprintf(fppts2,"%10d %10d\n",numbpts, numbtris);

for (i = 0; i < numbpts; i++) {
    fscanf(fppts1,"%lf %lf %lf %lf %lf %lf\n",
                                &x, &y, &z, &nx, &ny, &nz);
    
    det = 0.2*(z+5.0)*(z+5.0)/64.0;
    fprintf(fppts2, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
            x,y,z,x+det*nx, y+det*ny, z+det*nz, nx,ny,nz, nx,ny,nz);

}

for (i = 0; i < numbtris; i++) {
    fscanf(fppts1,"%d %d %d\n", &ii,&jj,&kk); 
    fprintf(fppts2, "%d %d %d\n", ii,jj,kk);
}

fclose(fppts1);
fclose(fppts2);
}


int CutTriangleByOnePlane(dir,p1,p2,p3,n1,n2,n3,a,tri1,tri2)
int   dir;
float p1[],p2[],p3[],n1[],n2[],n3[],a;
Triangle tri1, tri2;
{
int   i;
float q1[3], q2[3], q3[3], nn1[3],nn2[3], nn3[3],t12,t23, t13; 

if ( (p1[dir] >= a) && (p2[dir] >= a) && (p3[dir] >= a) ) {
   return(0);
}
if ( (p1[dir] <= a) && (p2[dir] <= a) && (p3[dir] <= a) ) {
   for (i = 0; i < 3; i++) {
      tri1->Vertex1[i] = p1[i];
      tri1->Vertex2[i] = p2[i];
      tri1->Vertex3[i] = p3[i];
      tri1->Normal1[i] = n1[i];
      tri1->Normal2[i] = n2[i];
      tri1->Normal3[i] = n3[i];
   }
   return(1);
}

if ( (p1[dir] <= p2[dir]) && (p2[dir] <= p3[dir]) ) {
  for (i = 0; i < 3; i++) {
      q1[i] = p1[i];
      q2[i] = p2[i];
      q3[i] = p3[i];
      nn1[i] = n1[i];
      nn2[i] = n2[i];
      nn3[i] = n3[i];
  }
}

if ( (p1[dir] <= p3[dir]) && (p3[dir] <= p2[dir]) ) {
  for (i = 0; i < 3; i++) {
      q1[i] = p1[i];
      q2[i] = p3[i];
      q3[i] = p2[i];
      nn1[i] = n1[i];
      nn2[i] = n3[i];
      nn3[i] = n2[i];
  }
}

if ( (p2[dir] <= p3[dir]) && (p3[dir] <= p1[dir]) ) {
  for (i = 0; i < 3; i++) {
      q1[i] = p2[i];
      q2[i] = p3[i];
      q3[i] = p1[i];
      nn1[i] = n2[i];
      nn2[i] = n3[i];
      nn3[i] = n1[i];
  }
}

if ( (p2[dir] <= p1[dir]) && (p1[dir] <= p3[dir]) ) {
  for (i = 0; i < 3; i++) {
      q1[i] = p2[i];
      q2[i] = p1[i];
      q3[i] = p3[i];
      nn1[i] = n2[i];
      nn2[i] = n1[i];
      nn3[i] = n3[i];
  }
}

if ( (p3[dir] <= p2[dir]) && (p2[dir] <= p1[dir]) ) {
  for (i = 0; i < 3; i++) {
      q1[i] = p3[i];
      q2[i] = p2[i];
      q3[i] = p1[i];
      nn1[i] = n3[i];
      nn2[i] = n2[i];
      nn3[i] = n1[i];
  }
}

if ( (p3[dir] <= p1[dir]) && (p1[dir] <= p2[dir]) ) {
  for (i = 0; i < 3; i++) {
      q1[i] = p3[i];
      q2[i] = p1[i];
      q3[i] = p2[i];
      nn1[i] = n3[i];
      nn2[i] = n1[i];
      nn3[i] = n2[i];
  }
}
if (q2[dir] >= a) {
   t12 = (a - q2[dir])/(q1[dir] - q2[dir]);
   t13 = (a - q3[dir])/(q1[dir] - q3[dir]);

   for (i = 0; i < 3; i++) {
       tri1->Vertex1[i] = q1[i];
       tri1->Normal1[i] = nn1[i];
       tri1->Vertex2[i] = t12*q1[i] + (1.0 - t12)*q2[i];
       tri1->Normal2[i] = t12*nn1[i] + (1.0 - t12)*nn2[i];
       tri1->Vertex3[i] = t13*q1[i] + (1.0 - t13)*q3[i];
       tri1->Normal3[i] = t13*nn1[i] + (1.0 - t13)*nn3[i];
   }
   return(1);
}
if (q2[dir] < a) {
   t23 = (a - q3[dir])/(q2[dir] - q3[dir]);
   t13 = (a - q3[dir])/(q1[dir] - q3[dir]);

   for (i = 0; i < 3; i++) {
       tri1->Vertex1[i] = q1[i];
       tri1->Normal1[i] = nn1[i];
       tri1->Vertex2[i] = q2[i];
       tri1->Normal2[i] = nn2[i];
       tri1->Vertex3[i] = t23*q2[i] + (1.0 - t23)*q3[i];
       tri1->Normal3[i] = t23*nn2[i] + (1.0 - t23)*nn3[i];

       tri2->Vertex1[i] = q1[i];
       tri2->Normal1[i] = nn1[i];
       tri2->Vertex2[i] = t23*q2[i] + (1.0 - t23)*q3[i]; 
       tri2->Normal2[i] = t23*nn2[i] + (1.0 - t23)*nn3[i];
       tri2->Vertex3[i] = t13*q1[i] + (1.0 - t13)*q3[i];
       tri2->Normal3[i] = t13*nn1[i] + (1.0 - t13)*nn3[i];

   }
   return(2);
}
return(0);
}
