#include <stdio.h>
#include <math.h>

#if ! defined(__APPLE__)
#include <malloc.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <Decimation/Linear/display.h>

static int  face[12][6] = {
   {1, 2, 3, 9, 4, 8},
   {2, 3, 0, 9, 5, 11},
   {3, 0, 1, 11,6, 10},
   {0, 1, 2, 10,7, 8},
   {5, 6, 7, 8, 0, 9},
   {6, 7, 4, 11,1, 9},
   {7, 4, 5, 10,2, 11},
   {4, 5, 6, 10,3, 8},
   {7, 10,3, 0, 9, 4},
   {5, 11,1, 4, 8, 0},
   {2, 11,6, 3, 8, 7},
   {1, 9, 5, 6, 10,2}
};

/* the edge index of the cube 

    c2/---------6-------/            
   / /:                /| b2
  /10 :              11 | |  
 / /  :              /  | |
c1/--------2--------/   | |
  |   7             |   5 |
  |   :             |   | |
  |   :             |   | b1
  3   :........4....1...|
  |  /              |  /
  | 8               | 9
  |/                |/
  /--------0--------/
  a1 ------------ a2
*/ 


/*---------------------------------------------------------------------------
PartialxyBernstein2D -- compute partial derivative about x,or y  of
                        a 2D Bernstein polynomial of degree m,n
----------------------------------------------------------------------------*/
void
PartialxyBernstein2D(weight,m,n,a1,a2,b1,b2,weight1,dir)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree m,n                              */
int      m,n;                    /* degrees in x and y directions           */
double a1,a2,b1,b2;              /* a1 <= x <= a2; b1 <= y <= b2;           */
double   weight1[];              /* the coefficients of derivatives about x,
                                    or y                                    */
int      dir;                    /* dir = 0--x, 1--y                        */
{
int      i,j,m1,n1;
double   x,y;

x = m/(a2 - a1);   y = n/(b2 - b1);
m1 = m + 1;  n1 = n + 1;

if (dir == 0) {
/* partial derivative about x */
for (i = 0; i < m; i++)
    for (j = 0; j < n1; j ++)
        weight1[n1*i + j] = x*(weight[n1*(i+1) + j] -
                               weight[n1*i + j]);
return;
}
if (dir == 1) {
/* partial derivative about y */
for (i = 0; i < m1; i++)
    for (j = 0; j < n; j ++)
        weight1[n*i + j] = y*(weight[n1*i + (j+1)] -
                              weight[n1*i + j]);
return;
}
}

/*---------------------------------------------------------------------------
PartialBernstein -- compute partial derivative about x,y and z of a 3D
                    Bernstein polynomial over a cube  of degree m,n,q
----------------------------------------------------------------------------*/
void
PartialBernstein(weight,m,n,q,a1,a2,b1,b2,c1,c2,dxweight,dyweight, dzweight)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree m,n,q                            */
int      m,n,q;                  /* the degree of the  Bernstein polynomial */
double a1,a2,b1,b2,c1,c2;        /* a1 <= x <= a2; b1 <= y <= b2;
                                    c1 <= z <= c2, the given cube           */
double   dxweight[],dyweight[], dzweight[];
                                 /* the coefficients of derivatives about x,
                                    y and z, respectively                   */
{
PartialxyzBernstein(weight,m,n,q,a1,a2,b1,b2,c1,c2,dxweight,0);
PartialxyzBernstein(weight,m,n,q,a1,a2,b1,b2,c1,c2,dyweight,1);
PartialxyzBernstein(weight,m,n,q,a1,a2,b1,b2,c1,c2,dzweight,2);
}

/*---------------------------------------------------------------------------
PartialxyzBernstein -- compute partial derivative about x,or y or z of  a 3D
                     Bernstein polynomial over a cube of degree m,n,q
----------------------------------------------------------------------------*/
void
PartialxyzBernstein(weight,m,n,q,a1,a2,b1,b2,c1,c2,weight1,dir)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree m,n,q                            */
int      m,n,q;                  /* the degree of the  Bernstein polynomial */
double a1,a2,b1,b2,c1,c2;        /* a1 <= x <= a2; b1 <= y <= b2;
                                    c1 <= z <= c2, the given cube           */
double   weight1[];
                                 /* the coefficients of derivatives about x,
                                    y or  z                                 */
int      dir;                    /* dir = 0--x, 1--y, 2--z                  */
{
int      i,j,k,m1,n1,q1;
double   x,y,z;

x = m/(a2 - a1);   y = n/(b2 - b1);   z = q/(c2 - c1);
m1 = m + 1;  n1 = n + 1;  q1 = q + 1;

if (dir == 0) {
/* partial derivative about x */
for (i = 0; i < m; i++)
    for (j = 0; j < n1; j ++)
        for (k = 0; k < q1; k ++)
            weight1[n1*q1*i + q1*j + k] = x*(weight[n1*q1*(i+1) + q1*j + k] -
                                            weight[n1*q1*i + q1*j + k]);
return;
}
if (dir == 1) {
/* partial derivative about y */
for (i = 0; i < m1; i++)
    for (j = 0; j < n; j ++)
        for (k = 0; k < q1; k ++)
            weight1[n*q1*i + q1*j + k] = y*(weight[n1*q1*i + q1*(j+1) + k] -
                                           weight[n1*q1*i + q1*j + k]);
return;
}
if (dir == 2) {
/* partial derivative about z */
for (i = 0; i < m1; i++)
    for (j = 0; j < n1; j ++)
        for (k = 0; k < q; k ++)
            weight1[n1*q*i + q*j + k] = z*(weight[n1*q1*i + q1*j + k+1] -
                                            weight[n1*q1*i + q1*j + k]);
}
}
/*---------------------------------------------------------------------------
Intersection12Edge -- the intersection of 12 edges with a Bernsein
----------------------------------------------------------------------------*/
void
Intersection12Edge(weight,m,a1,a2,b1,b2,c1,c2,intersect,numberin,flag)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree (2,2,2) or  (3,3,3)              */
int      m;                      /* degree m = 2 or m = 3                   */
double a1,a2,b1,b2,c1,c2;        /* a1 <= x <= a2; b1 <= y <= b2;
                                    c1 <= z <= c2                           */
double intersect[12][2];         /* the intersection points                 */
int    *numberin;                /* the number of intersection points       */
int    *flag;                    /* condition a,b,c in paper are fail       */
{
double intersect1[12][2];
int    i,m1,m2;
double   weight1[64], max;

/*
printf("in Intersection12Edge\n");
print_vector(weight,8);
*/
for (i = 0; i < 12; i++){
    intersect[i][0] = 0.0;
    intersect[i][1] = 0.0;
    intersect1[i][0] = 0.0;
    intersect1[i][1] = 0.0;
}
m1 = m + 1;   m2 = m1*m1;
max = -10.0;
for (i = 0; i < m1*m2; i ++) {
    if (fabs(weight[i]) > max) max = fabs(weight[i]);
}
 
/*
for (i = 0; i < m1*m2; i ++) weight1[i] = weight[i] + 0.000314159261828;
*/
for (i = 0; i < m1*m2; i ++) weight1[i] = weight[i];

Intersection12Edge1(weight1,m,a1,a2,b1,b2,c1,c2,intersect,numberin,flag);
if (PointsNum(intersect)>6  || *flag == 0) {
   *flag = 0;
   return;
}

/*
Intersection12Edge1(weight,m,a1,a2,b1,b2,c1,c2,intersect1,numberin,flag);
*flag = 1;
for (i = 0; i < 12; i ++) intersect[i][1] = intersect1[i][1];
*/
*numberin = PointsNum(intersect);

if (intersect[0][0]+intersect[8][0]+intersect[4][0]+intersect[9][0] > 2.0){
   *flag = 0;
   return;
}
if (intersect[2][0]+intersect[6][0]+intersect[10][0]+intersect[11][0] > 2.0){
   *flag = 0;
   return;
}
if (intersect[0][0]+intersect[1][0]+intersect[2][0]+intersect[3][0] > 2.0){
   *flag = 0;
   return;
}
if (intersect[4][0]+intersect[5][0]+intersect[6][0]+intersect[7][0] > 2.0){
   *flag = 0;
   return;
}
if (intersect[1][0]+intersect[5][0]+intersect[9][0]+intersect[11][0] > 2.0){
   *flag = 0;
   return;
}
if (intersect[3][0]+intersect[7][0]+intersect[8][0]+intersect[10][0] > 2.0){
   *flag = 0;
   return;
}

}

/*---------------------------------------------------------------------------
Intersection12Edge1 -- the intersection of 12 edges with a Bernsein
----------------------------------------------------------------------------*/
void
Intersection12Edge1(weight,m,a1,a2,b1,b2,c1,c2,intersect,numberin,flag)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree (2,2,2) or  (3,3,3)              */
int      m;                      /* degree m = 2 or m = 3                   */
double a1,a2,b1,b2,c1,c2;        /* a1 <= x <= a2; b1 <= y <= b2;
                                    c1 <= z <= c2                           */
double intersect[12][2];         /* the intersection points                 */
int    *numberin;                /* the number of intersection points       */
int    *flag;                    /* condition a,b,c in paper are fail       */
{
double edges[12][4];
int    i,j,k,sum,m1,m2,flag1;

m1 = m + 1;   m2 = m1*m1;
for (i = 0; i < m1; i ++){
    edges[0][i] = weight[m2*i];
    edges[1][i] = weight[m*m2 +m1*i];
    edges[2][i] = weight[m2*i + m*m1];
    edges[3][i] = weight[m1*i];
    edges[4][i] = weight[m2*i+m];
    edges[5][i] = weight[m*m2+m1*i+m];
    edges[6][i] = weight[m2*i+ m*m1+m];
    edges[7][i] = weight[m1*i+m];
    edges[8][i] = weight[i];
    edges[9][i] = weight[m*m2+i];
    edges[10][i] = weight[m*m1+i];
    edges[11][i] = weight[m*m2+m*m1+i];
}
sum = 0;
*flag = 1;
for (i = 0; i < 12; i ++){
    FindRoots23Edges(edges,i,intersect,&flag1,m);
    sum = sum  + intersect[i][0];
    if (flag1 == 0) *flag = 0;
}
*numberin = sum;
}

/*---------------------------------------------------------------------------
Intersection4Edge -- the intersection of 12 edges with a Bernsein
----------------------------------------------------------------------------*/
void
Intersection4Edge(weight,m,a1,a2,b1,b2,intersect,numberin,flag)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree (2,2,2) or  (3,3,3)              */
int      m;                      /* degree m = 2 or m = 3                   */
double a1,a2,b1,b2;              /* a1 <= x <= a2; b1 <= y <= b2;           */
double intersect[4][2];          /* the intersection points                 */
int    *numberin;                /* the number of intersection points       */
int    *flag;                    /* condition a,b,c in paper are fail       */
{
double intersect1[4][2];
int    i,sum,m1;
double   weight1[16];

for (i = 0; i < 4; i++){
    intersect[i][0] = 0.0;
    intersect[i][1] = 0.0;
    intersect1[i][0] = 0.0;
    intersect1[i][1] = 0.0;
}

m1 = m + 1;
sum = 0.0;
for (i = 0; i < 16; i ++) weight1[i] = weight[i];
sum = 0;
if (weight[0] == 0.0) {
   weight1[0] = 0.000001;
   sum = 1;
}
if (weight[m1*m] == 0.0) {
   weight1[m1*m] = 0.000001;
   sum = 1;
}
if (weight[m] == 0.0) {
   weight1[m] = 0.000001;
   sum = 1;
}
if (weight[m+m1*m] == 0.0) {
   weight1[m +m1*m] = 0.000001;
   sum = 1;
}
Intersection4Edge1(weight1,m,a1,a2,b1,b2,intersect,numberin,flag);
if (PointsNum2D(intersect) > 2  || *flag == 0) {
   *flag = 0;
   return;
}

if (sum > 0){
  Intersection4Edge1(weight,m,a1,a2,b1,b2,intersect1,numberin,flag);
  *flag = 1;
  for (i = 0; i < 4; i ++) intersect[i][1] = intersect1[i][1];
  *numberin = PointsNum2D(intersect);
}

}
/*---------------------------------------------------------------------------
PointsNum2D -- the number of the zero points
----------------------------------------------------------------------------*/
int
PointsNum2D(intersect)
double  intersect[4][2];
{
int i,j;
j = 0;
for (i = 0; i < 4; i++){
    if (intersect[i][0] == 1.0) j = j + 1;
}
return j;
}
/*---------------------------------------------------------------------------
Intersection4Edge1 -- the intersection of 12 edges with a Bernsein
----------------------------------------------------------------------------*/
void
Intersection4Edge1(weight,m,a1,a2,b1,b2,intersect,numberin,flag)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree (2,2,2) or  (3,3,3)              */
int      m;                      /* degree m = 2 or m = 3                   */
double a1,a2,b1,b2;              /* a1 <= x <= a2; b1 <= y <= b2;           */
double intersect[4][2];          /* the intersection points                 */
int    *numberin;                /* the number of intersection points       */
int    *flag;                    /* condition a,b,c in paper are fail       */
{
double edges[4][4];
int    i,j,k,sum,m1,flag1;

/*
print_vector(weight,16);
*/
m1 = m + 1;
for (i = 0; i < m1; i ++){
    edges[0][i] = weight[m1*i];
    edges[1][i] = weight[m*m1 +i];
    edges[2][i] = weight[m1*i + m];
    edges[3][i] = weight[i];
}
sum = 0;
*flag = 1;
for (i = 0; i < 4; i ++){
    FindRoots2D123Edges(edges,i,intersect,&flag1,m);
    sum = sum  + intersect[i][0];
    if (flag1 == 0) {
       *flag = 0;
    }
}
*numberin = sum;
}
/*---------------------------------------------------------------------------
FindRoots23Edges -- find roots of 12 edges of degree 2 or 3
----------------------------------------------------------------------------*/
void
FindRoots23Edges(edges,i,intersect,flag,degree)
double edges[12][4];
int    i, *flag;
double intersect[12][2];
int    degree;
{
double  weight[4],roots[3];
int     j,k;

for (j = 0; j < degree + 1; j ++)
    weight[j] = edges[i][j];
BernsteinRoot(weight,roots,&k,degree);
if (k > 1) {
   *flag = 0;
   intersect[i][0] = k;
   return;
}
if (k == 1) {
   *flag = 1;
   intersect[i][0] = k;
   intersect[i][1] = roots[0];
   return;
}
if (k == 0) {
   *flag = 1;
   intersect[i][0] = k;
}
}
/*---------------------------------------------------------------------------
FindRoots2D123Edges -- find roots of 4 edges of degree 1 or 2 or 3
----------------------------------------------------------------------------*/
void
FindRoots2D123Edges(edges,i,intersect,flag,degree)
double edges[4][4];
int    i, *flag;
double intersect[4][2];
int    degree;
{
double  weight[4],roots[3];
int     j,k;

for (j = 0; j < degree + 1; j ++)
    weight[j] = edges[i][j];
BernsteinRoot(weight,roots,&k,degree);
if (k > 1) {
   *flag = 0;
   intersect[i][0] = k;
   return;
}
if (k == 1) {
   *flag = 1;
   intersect[i][0] = k;
   intersect[i][1] = roots[0];
   return;
}
if (k == 0) {
   *flag = 1;
   intersect[i][0] = k;
}
}

/*---------------------------------------------------------------------------
DisplayBerns123 -- linearize a iso-surface befined by c = f(x,y,z), 
                   here f is a Bernstein polynomial of degree 
                   (1,1,1) or (2,2,2) or (3,3,3) over a cube 
                   [a1,a2]*[b1,b2]*[c1,c2] by a recursively 
                   subdivision, where c is a given threshold value
----------------------------------------------------------------------------*/
void
DisplayBerns123(weight,a1,a2,b1,b2,c1,c2, m,epsilon,scaler,threshold,polygons)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree (1,1,1) or (2,2,2) or  (3,3,3)   */
double a1,a2,b1,b2,c1,c2;        /* a1 <= x <= a2; b1 <= y <= b2;
                                    c1 <= z <= c2, the given cube           */
int    m;                        /* degree,  m = 1 or m = 2, or m = 3       */
double epsilon;                  /* the upper bound of the partition size   */ 
double scaler[];                 /* the mapping parameter                   */
double threshold;                /* the threshold value                     */
Polygons *polygons;              /* polygons chain  created                 */
{
double   dxweight[48],dyweight[48], dzweight[48];
double   cube[6],weight2[64];
int      i,m1;
Polygons endpoly;

if (*polygons != NULL) {
   fprintf(stderr,
        "polygons need to be NUll pointer, now it is force to be NULL!\n");
   *polygons = NULL;
}
m1 = (m+1)*(m+1)*(m+1);
for (i = 0; i < m1; i++) weight2[i] = weight[i] - threshold + 0.000314159261828;;
cube[0] = a1;    cube[1] = a2;   cube[2] = b1;
cube[3] = b2;    cube[4] = c1;   cube[5] = c2;
PartialBernstein(weight2,m,m,m,a1,a2,b1,b2,c1,c2,dxweight,dyweight, dzweight);
DisplayBerns123in(weight2,dxweight,dyweight, dzweight,a1,a2,b1,b2,c1,c2,
                  m,epsilon,cube,scaler,polygons,&endpoly);
}

/*---------------------------------------------------------------------------
DisplayBernsL123 -- linearize a iso-surface befined by c = f(x,y,z), 
                   here f is a Bernstein polynomial of degree 
                   (1,1,1) or (2,2,2) or (3,3,3) over a cube 
                   [a1,a2]*[b1,b2]*[c1,c2] by a recursively 
                   subdivision, where c is a given threshold value
----------------------------------------------------------------------------*/
void
DisplayBernsL123(weight,a1,a2,b1,b2,c1,c2, m,epsilon,scaler,threshold,polygons)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree (1,1,1) or (2,2,2) or  (3,3,3)   */
double a1,a2,b1,b2,c1,c2;        /* a1 <= x <= a2; b1 <= y <= b2;
                                    c1 <= z <= c2, the given cube           */
int    m;                        /* degree,  m = 1 or m = 2, or m = 3       */
double epsilon;                  /* the upper bound of the partition size   */ 
double scaler[];                 /* the mapping parameter                   */
double threshold;                /* the threshold value                     */
PolyLink polygons;               /* polygons chain  created                 */
{
double   dxweight[48],dyweight[48], dzweight[48];
double   cube[6],weight2[64];
int      i,m1;

m1 = (m+1)*(m+1)*(m+1);
for (i = 0; i < m1; i++) weight2[i] = weight[i] - threshold;
cube[0] = a1;    cube[1] = a2;   cube[2] = b1;
cube[3] = b2;    cube[4] = c1;   cube[5] = c2;
PartialBernstein(weight2,m,m,m,a1,a2,b1,b2,c1,c2,dxweight,dyweight, dzweight);
DisplayBerns123in(weight2,dxweight,dyweight, dzweight,a1,a2,b1,b2,c1,c2,
                  m,epsilon,cube,scaler,&polygons->head,&polygons->tail);
}

/*---------------------------------------------------------------------------
DisplayBerns123in -- display a implicitly defined surface befined by 
                     0 = f(x,y,z), here f is a Bernstein polynomial 
                     of degree (1,1,1) or (2,2,2) or (3,3,3) over a
                     cube [a1,a2]*[b1,b2]*[c1,c2] by a recursively 
                     subdivision
----------------------------------------------------------------------------*/
void
DisplayBerns123in(weight,dxweight,dyweight, dzweight,a1,a2,b1,b2,c1,c2,
				m,epsilon,cube,scaler,polygons,endpoly)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree (1,1,1) or (2,2,2) or  (3,3,3)   */
double   dxweight[],dyweight[], dzweight[];
                                 /* the coefficients of derivatives about x,
                                    y and z                                 */
double a1,a2,b1,b2,c1,c2;        /* a1 <= x <= a2; b1 <= y <= b2;
                                    c1 <= z <= c2                           */
int    m;                        /* degree m = 2, or m = 3                  */
double epsilon;                  /* the size conttrol of the partition      */
double cube[],scaler[];
Polygons *polygons,*endpoly;
{
int    i,j,k,sum,sum1,m1,m2,flag,numberin,dir;
double aa1,aa2,bb1,bb2,cc1,cc2,a1a1,a1a2,b1b1,b1b2,c1c1,c1c2;
double   weight1[64], weight2[64];
double intersect[12][2];         /* the intersection points                 */

m1 = m + 1;
m2 = m1*m1*m1;

/*
printf("intervals %f,%f,%f,%f,%f,%f\n", a1,a2,b1,b2,c1,c2);
print_vector(weight,8);
*/

/* no surface test         */
sum = 0;  sum1 = 0;
for (i = 0; i < m2; i++) {
    if (weight[i] >= 0.0 ) sum = sum + 1;
    if (weight[i] <= 0.0 ) sum1 = sum1 - 1;
}
if (sum  == m2 || sum1 == -m2) {
   return;
}
/* partition direction    */
   aa1 = a1;  aa2 = a2;   bb1 = b1;  bb2 = b2;  cc1 = c1;  cc2 = c2;
   a1a1 = a1; a1a2 = a2;  b1b1 = b1; b1b2 = b2; c1c1 = c1; c1c2 = c2;
if ((a2 - a1 >= b2 - b1) && (a2 - a1 >= c2 - c1) ) {
   dir = 1;
   aa2 = (a1 + a2)/2.0;
   a1a1 = aa2;
} else
if ((b2 - b1 >= a2 - a1) && (b2 - b1 >= c2 - c1) ){
   dir = 2;
   bb2 = (b1 + b2)/2.0;
   b1b1 = bb2;
} else
if ((c2 - c1 >= a2 - a1) && (c2 - c1 >= b2 - b1) ) {
   dir = 3;
   cc2 = (c1 + c2)/2.0;
   c1c1 = cc2;
}

/* partition              */
if ((a2 - a1 > epsilon) || (b2 - b1 >epsilon) || (c2 - c1 > epsilon) ){
   Partition(weight,dir,m,0.5,weight1,weight2);
   DisplayBerns123in(weight1,dxweight,dyweight, dzweight,aa1,aa2,bb1,bb2,
                     cc1,cc2,m,epsilon,cube,scaler,polygons,endpoly);
   DisplayBerns123in(weight2,dxweight,dyweight, dzweight,a1a1,a1a2,b1b1,b1b2,
                     c1c1,c1c2,m,epsilon,cube,scaler,polygons,endpoly);
}  else {
   Intersection12Edge(weight,m,a1,a2,b1,b2,c1,c2,intersect,&numberin,&flag);
/*
printf("&numberin,&flag %d,%d\n",numberin,flag);
*/
   if (flag == 0) {
      Partition(weight,dir,m,0.5,weight1,weight2);
      DisplayBerns123in(weight1,dxweight,dyweight, dzweight,aa1,aa2,bb1,bb2,
                        cc1,cc2,m,epsilon,cube,scaler,polygons,endpoly);
      DisplayBerns123in(weight2,dxweight,dyweight, dzweight, a1a1,a1a2,b1b1,b1b2,
                        c1c1,c1c2,m,epsilon,cube,scaler,polygons,endpoly);
   }
   if (flag == 1 && numberin >= 3) {
      DisplayPolygons(intersect,dxweight,dyweight, dzweight,m,
                      a1,a2,b1,b2,c1,c2,cube,scaler,polygons,endpoly);
   }
}

}

/*---------------------------------------------------------------------------
DisplayPolygons -- display polygons
----------------------------------------------------------------------------*/
void
DisplayPolygons(intersect,dxweight,dyweight, dzweight,m,
                 a1,a2,b1,b2,c1,c2,cube,scaler,polygons,endpoly)
double intersect[12][2];
double dxweight[],dyweight[], dzweight[];
int    m;
double a1,a2,b1,b2,c1,c2;        /* a1 <= x <= a2; b1 <= y <= b2;
                                    c1 <= z <= c2                           */
double cube[6],scaler[];
Polygons *polygons, *endpoly;
{
double x,y,z,point[3],point1[3];
int    i,start,j,k,sum,run,new,old;
Vertex  *vertices;
Normal  *normals;
double intersect1[12][2];
Polygons temp;

/*
printf("PointsNum %d \n", PointsNum(intersect));
print_vector(intersect,24);

print_vector(cube,6);
printf("%f,%f,%f,%f,%f,%f\n",a1,a2,b1,b2,c1,c2);
*/
for (i = 0; i < 12; i++) {
    intersect1[i][0] = intersect[i][0];
    intersect1[i][1] = intersect[i][1];
}

/* temarary change for -----*/
for (i = 0; i < 12; i++) {
    if (intersect[i][0] == 1.0){
       if (intersect[i][1] < 0.0000001 || intersect[i][1] > 0.9999999) {
       /*
          intersect[i][0] = 0.0;
          */
          printf("what happen = %f\n", intersect[i][1]);
          print_vector(intersect1,24);
          /*
          return;
          */
       }
    }
}
/* end temarary change */
if (PointsNum(intersect) < 3) return;

vertices = (Vertex *) malloc (6 * sizeof (Vertex));
normals = (Normal *) malloc (6 * sizeof (Normal));

run = -1; start = -1;   sum = 0;

for (i = 0; i < 12 && start < 0; i++){
    if (intersect[i][0] == 1.0) start = i;
}
new = start;  old = -1;
while(start != run) {
    if (PointsNum(intersect) == 0.0 || sum > 5) {
       printf("---------------------------------------------------------\n");
       print_vector(intersect1,24);
       return;
    }
    if (which(new,old) < 3 ){
         for (i = 0; i < 3; i++){
            printf("face = %d\n", face[new][i + 3]);
            if (intersect[face[new][i + 3]][0] == 1.0 ){
                run = face[new][i + 3];
                intersect[face[new][i + 3]][0] = 0.0;
            }
         }
    } else
    {
         for (i = 0; i < 3; i++){
            printf("face = %d\n", face[new][i]);
            if (intersect[face[new][i]][0] == 1.0 ){
                run = face[new][i];
                intersect[face[new][i]][0] = 0.0;
            }
         }
    }
    old = new;        new = run;
/*
    printf("new,old, sum %d,%d,%d\n", new,old,sum);
*/
    GotPoint(intersect,a1,a2,b1,b2,c1,c2,new,point,point1,scaler);
    for (i = 0; i < 3 ; i++) vertices[sum][i] = point[i];

    x = (point1[0] - cube[0])/(cube[1] - cube[0]);
    y = (point1[1] - cube[2])/(cube[3] - cube[2]);
    z = (point1[2] - cube[4])/(cube[5] - cube[4]);
    normals[sum][0] = BernsteinValue(dxweight,x,y,z,m-1,m,m);
    normals[sum][1] = BernsteinValue(dyweight,x,y,z,m,m-1,m);
    normals[sum][2] = BernsteinValue(dzweight,x,y,z,m,m,m-1);

    x = 0.0;
    for (i = 0; i < 3 ; i++) x = x + normals[sum][i]*normals[sum][i];
    x = sqrt(x);
    for (i = 0; i < 3 ; i++) normals[sum][i] = normals[sum][i]/x;
    sum = sum + 1;
}
if (cube[3] - cube[2] < 0.0002 && cube[5] - cube[4] < 0.0002) {
   /*
   LabelledPointAbs(vertices[0][0], vertices[0][1],vertices[0][2]);
   */ 
   temp = CreatPolygons(1, vertices, NULL);
   if (*polygons == NULL){
      *polygons = temp;
      *endpoly = temp;
   }  else {
      (*endpoly)->next = temp;
      *endpoly = temp;
   }
}  else {
   temp = CreatPolygons(sum, vertices, normals);
   if (*polygons == NULL){
      *polygons = temp;
      *endpoly = temp;
   }  else {
      (*endpoly)->next = temp;
      *endpoly = temp;
   }
   /*
   DrawPolygon(sum, vertices, normals);
   */
}
if (PointsNum(intersect) >= 3 )
   DisplayPolygons(intersect,dxweight,dyweight, dzweight,m,
                 a1,a2,b1,b2,c1,c2,cube,scaler,polygons,endpoly);
}
/*---------------------------------------------------------------------------
DisplayPolygons2D -- display polygons
----------------------------------------------------------------------------*/
void
DisplayPolygons2D(intersect,m, a1,a2,b1,b2,cube,zvalue,scaler,
                  polygons,endpoly)
double intersect[4][2];
int    m;
double a1,a2,b1,b2;             /* a1 <= x <= a2; b1 <= y <= b2;            */
double cube[6],zvalue,scaler[];
Polygons  *polygons, *endpoly;
{
double x,y,z,x1,y1,z1;
int    i,start,end;
Vertex  *vertices;
Polygons  temp;
/* temarary change for -----*/
/*
for (i = 0; i < 4; i++) {
    if (intersect[i][0] == 1.0){
       if (intersect[i][1] < 0.0000001 || intersect[i][1] > 0.9999999) {
          intersect[i][0] = 0.0;
          printf("what happen = %f\n", intersect[i][1]);
          return;
       }
    }
}
*/
/* end temarary change */

start = -1;
for (i = 0; i < 4; i++){
    if (intersect[i][0] == 1.0 && start == -1) start = i;
    if (intersect[i][0] == 1.0 && start > -1 ) end  = i;
}
if (start == 0) { x = intersect[start][1]*(a2-a1) + a1; y = b1;   z = zvalue; }
if (start == 1) { x = a2; y = intersect[start][1]*(b2-b1) + b1;   z = zvalue; }
if (start == 2) { x = intersect[start][1]*(a2-a1) + a1; y = b2;   z = zvalue; }
if (start == 3) { x = a1; y = intersect[start][1]*(b2-b1) + b1;   z = zvalue; }
x = scaler[0]*x + scaler[1];
y = scaler[0]*y + scaler[2];
z = scaler[0]*z + scaler[3];


if (end == 0) { x1 = intersect[end][1]*(a2-a1) + a1; y1 = b1;   z1 = zvalue; }
if (end == 1) { x1 = a2; y1 = intersect[end][1]*(b2-b1) + b1;   z1 = zvalue; }
if (end == 2) { x1 = intersect[end][1]*(a2-a1) + a1; y1 = b2;   z1 = zvalue; }
if (end == 3) { x1 = a1; y1 = intersect[end][1]*(b2-b1) + b1;   z1 = zvalue; }
x1 = scaler[0]*x1 + scaler[1];
y1 = scaler[0]*y1 + scaler[2];
z1 = scaler[0]*z1 + scaler[3];
vertices = (Vertex *) malloc (2 * sizeof (Vertex));
vertices[0][0] = x;    vertices[0][1] = y;    vertices[0][2] = z;
vertices[1][0] = x1;   vertices[1][1] = y1;   vertices[1][2] = z1;
/*
LineAbs(x,y,z,x1,y1,z1);
*/
temp = CreatPolygons(2, vertices, NULL);
if (*polygons == NULL){
   *polygons = temp;
   *endpoly = temp;
}  else {
   (*endpoly)->next = temp;
   *endpoly = temp;
}
}
/*---------------------------------------------------------------------------
DisplayBerns2D123 -- linearize  a iso-curve befined by c = f(x,y),
                   at the plane defined z = constant, 
                   here f is a Bernstein polynomial of degree
                   (1,1) or (2,2) or (3,3) over a rectangular 
                   [a1,a2]*[b1,b2] by a recursively subdivision,
                   where c is a given threshold value
----------------------------------------------------------------------------*/
void
DisplayBerns2D123(weight,a1,a2,b1,b2,m,epsilon,zz,scaler,threshold,polygons)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree (1,1) or (2,2) or  (3,3)         */
double a1,a2,b1,b2;              /* a1 <= x <= a2; b1 <= y <= b2; the given
                                    rectangular                             */
int    m;                        /* degree,   m = 1 or m = 2, or m = 3      */
double epsilon;                  /* the size conttrol of the partition      */
double zz;                       /* the z-value of the plane z = zz on which
                                    the iso-curve will be displayed.        */
double scaler[];                 /* the mapping parameter                   */
double threshold;                /* the threshold value                     */
Polygons *polygons;              /* polygons created                        */
{
double   cube[6],weight2[16];
int      i,m1;
Polygons endpoly;

if (*polygons != NULL) {
   fprintf(stderr,
        "polygons need to be NUll pointer, now it is force to be NULL!\n");
   *polygons = NULL;
}
m1 = (m+1)*(m+1);
for (i = 0; i < m1; i++) weight2[i] = weight[i] - threshold;
cube[0] = a1;    cube[1] = a2;   cube[2] = b1;
cube[3] = b2;   
DisplayBerns2D123in(weight2,a1,a2,b1,b2,m,epsilon,cube,zz,scaler,
                    polygons,&endpoly);
}


/*---------------------------------------------------------------------------
DisplayBerns2D123in -- display 2D Bernstein polynomial degree (1,1) or
                     (2,2) or (3,3) recursively
----------------------------------------------------------------------------*/
void
DisplayBerns2D123in(weight,a1,a2,b1,b2,m,epsilon,cube,zvalue,scaler,
                    polygons,endpoly)
double   weight[];               /* the coefficients of Bernstein polynomial
                                    degree  (1,1) or  (2,2) or  (3,3)       */
double a1,a2,b1,b2;              /* a1 <= x <= a2; b1 <= y <= b2;           */
int    m;                        /* degree m = 1 or m = 2, or m = 3         */
double epsilon;                  /* the size conttrol of the partition      */
double cube[],zvalue,scaler[];
Polygons *polygons,*endpoly;
{
int    i,j,k,sum,sum1,m1,m2,flag,numberin,dir;
double aa1,aa2,bb1,bb2,a1a1,a1a2,b1b1,b1b2;
double   weight1[16], weight2[16];
double intersect[4][2];          /* the intersection points                 */

m1 = m + 1;
m2 = m1*m1;

/* no curve  test         */
sum = 0;  sum1 = 0;

for (i = 0; i < m2; i++) {
    if (weight[i] >= 0.0 ) sum = sum + 1;
    if (weight[i] <= 0.0 ) sum1 = sum1 - 1;
}
if (sum  == m2 || sum1 == -m2) return;

/* partition direction    */
   aa1 = a1;  aa2 = a2;   bb1 = b1;  bb2 = b2;
   a1a1 = a1; a1a2 = a2;  b1b1 = b1; b1b2 = b2;
if ((a2 - a1 >= b2 - b1) ) {
   dir = 1;
   aa2 = (a1 + a2)/2.0;
   a1a1 = aa2;
} else
if ((b2 - b1 >= a2 - a1) ){
   dir = 2;
   bb2 = (b1 + b2)/2.0;
   b1b1 = bb2;
}
/* partition              */
if ((a2 - a1 > epsilon) || (b2 - b1 >epsilon)  ){
   Partition2D(weight,dir,m,0.5,weight1,weight2);
   DisplayBerns2D123in(weight1, aa1,aa2,bb1,bb2,m,epsilon,cube,zvalue,
                       scaler,polygons,endpoly);
   DisplayBerns2D123in(weight2, a1a1,a1a2,b1b1,b1b2,m,epsilon,cube,zvalue,
                       scaler,polygons,endpoly);
}  else {
   Intersection4Edge(weight,m,a1,a2,b1,b2,intersect,&numberin,&flag);
/*
print_vector(weight,4);
print_vector(intersect,8);
printf("&numberin,&flag %d,%d\n",numberin,flag);
*/
   if (flag == 0) {
      Partition2D(weight,dir,m,0.5,weight1,weight2);
      DisplayBerns2D123in(weight1, aa1,aa2,bb1,bb2,m,epsilon,cube,zvalue,
                          scaler,polygons,endpoly);
      DisplayBerns2D123in(weight2, a1a1,a1a2,b1b1,b1b2,m,epsilon,cube,zvalue,
                          scaler,polygons,endpoly);
   }
   if (flag == 1 && numberin == 2) {
      DisplayPolygons2D(intersect,m, a1,a2,b1,b2,cube,zvalue,scaler,
                        polygons,endpoly);
   }
}

}

/*---------------------------------------------------------------------------
----------------------------------------------------------------------------*/
int
which(i,j)
int  i,j;
{
int  k,iw;
iw = -1;
for (k = 0; k < 6; k++){
    if (face[i][k] == j) iw = k;
}
return iw;
}

/*---------------------------------------------------------------------------
GotPoint -- get a point
----------------------------------------------------------------------------*/
void
GotPoint(intersect,a1,a2,b1,b2,c1,c2,i,point,point1,scaler)
double intersect[12][2];
double a1,a2,b1,b2,c1,c2;
int    i;
double point[],point1[],scaler[];
{
double x,y,z;

       if (i == 0) { x = intersect[i][1]*(a2-a1) + a1; y = b1;   z = c1; }
       if (i == 1) { x = a2; y = intersect[i][1]*(b2-b1) + b1;   z = c1; }
       if (i == 2) { x = intersect[i][1]*(a2-a1) + a1; y = b2;   z = c1; }
       if (i == 3) { x = a1; y = intersect[i][1]*(b2-b1) + b1;   z = c1; }

       if (i == 4) { x = intersect[i][1]*(a2-a1) + a1; y = b1;   z = c2; }
       if (i == 5) { x = a2; y = intersect[i][1]*(b2-b1) + b1;   z = c2; }
       if (i == 6) { x = intersect[i][1]*(a2-a1) + a1; y = b2;   z = c2; }
       if (i == 7) { x = a1; y = intersect[i][1]*(b2-b1) + b1;   z = c2; }

       if (i == 8) { x = a1; y = b1;   z = intersect[i][1]*(c2-c1) + c1; }
       if (i == 9) { x = a2; y = b1;   z = intersect[i][1]*(c2-c1) + c1; }
       if (i == 10) { x = a1; y = b2;   z = intersect[i][1]*(c2-c1) + c1; }
       if (i == 11) { x = a2; y = b2;   z = intersect[i][1]*(c2-c1) + c1; }
/*  scale                */
/*
printf("i,number, x,y,z = %d,%d %f,%f,%f\n",i,PointsNum(intersect), x,y,z);
*/
       point1[0] = x; point1[1] = y;  point1[2] = z;
       x = scaler[0]*x + scaler[1];
       y = scaler[0]*y + scaler[2];
       z = scaler[0]*z + scaler[3];

/*
       LabelledPointAbs(x,y,z);
*/
       point[0] = x; point[1] = y;  point[2] = z;

}

/*---------------------------------------------------------------------------
PointsNum -- the number of the zero points
----------------------------------------------------------------------------*/
int
PointsNum(intersect)
double  intersect[12][2];
{
int i,j;
j = 0;
for (i = 0; i < 12; i++){
    if (intersect[i][0] == 1.0) j = j + 1;
}
return j;
}

/*---------------------------------------------------------------------------
Partition -- partition a Berbstein Polynomial into two
----------------------------------------------------------------------------*/
void
Partition(weight,dir,m,c,weight1,weight2)
double   weight[], weight1[],weight2[];
int      dir,m;
double   c;
{
int      i,j,k,l,m1,m2;
double   sum;

m1 = m + 1;
m2 = m1*m1;
/* x-dircetion  */
for (i = 0; i < m1; i ++)
    for (j = 0; j < m1; j ++)
        for (k = 0; k < m1; k ++) {
           if ( dir == 1) {
               sum = 0.0;
               for (l = 0; l <= i; l ++){
                   sum = sum + weight[m2*l + j*m1 + k]*BernsteinBaseValue(i,l,c);
               }
               weight1[m2*i + j*m1 + k] = sum;
               sum = 0.0;
               for (l = 0; l <= m - i; l ++){
                   sum = sum + weight[m2*(i+l) + j*m1 + k]
                         *BernsteinBaseValue(m-i,l,c);
               }
               weight2[m2*i + j*m1 + k] = sum;
           }
           if ( dir == 2) {
              sum = 0.0;
               for (l = 0; l <= j; l ++){
                   sum = sum + weight[m2*i + l*m1 + k]*BernsteinBaseValue(j,l,c);
               }
               weight1[m2*i + j*m1 + k] = sum;
               sum = 0.0;
               for (l = 0; l <= m - j; l ++){
                   sum = sum + weight[m2*i + (l+j)*m1 + k]
                         *BernsteinBaseValue(m-j,l,c);
               }
               weight2[m2*i + j*m1 + k] = sum;
           }
           if ( dir == 3) {
              sum = 0.0;
               for (l = 0; l <= k; l ++){
                   sum = sum + weight[m2*i + j*m1 + l]*BernsteinBaseValue(k,l,c);
               }
               weight1[m2*i + j*m1 + k] = sum;
               sum = 0.0;
               for (l = 0; l <= m - k; l ++){
                   sum = sum + weight[m2*i + j*m1 + k+l]
                         *BernsteinBaseValue(m-k,l,c);
               }
               weight2[m2*i + j*m1 + k] = sum;
           }
       }
}

/*---------------------------------------------------------------------------
Partition2D -- partition a Berbstein Polynomial into two
----------------------------------------------------------------------------*/
void
Partition2D(weight,dir,m,c,weight1,weight2)
double   weight[], weight1[],weight2[];
int      dir,m;
double   c;
{
int      i,j,l,m1;
double   sum;

m1 = m + 1;
/* x-dircetion  */
for (i = 0; i < m1; i ++)
    for (j = 0; j < m1; j ++) {
           if ( dir == 1) {
               sum = 0.0;
               for (l = 0; l <= i; l ++){
                   sum = sum + weight[m1*l + j]*BernsteinBaseValue(i,l,c);
               }
               weight1[m1*i + j] = sum;
               sum = 0.0;
               for (l = 0; l <= m - i; l ++){
                   sum = sum + weight[m1*(i+l) + j]
                         *BernsteinBaseValue(m-i,l,c);
               }
               weight2[m1*i + j] = sum;
           }
           if ( dir == 2) {
              sum = 0.0;
               for (l = 0; l <= j; l ++){
                   sum = sum + weight[m1*i + l]*BernsteinBaseValue(j,l,c);
               }
               weight1[m1*i + j] = sum;
               sum = 0.0;
               for (l = 0; l <= m - j; l ++){
                   sum = sum + weight[m1*i + (l+j)]
                         *BernsteinBaseValue(m-j,l,c);
               }
               weight2[m1*i + j] = sum;
           }
   }
}

/*---------------------------------------------------------------------------
Display2DFunction --linearize the surface z = f(x,y), here f is in BB-form
                    of degree m,n over a rectangular [a1,a2]*[b1,b2].
                    This function can also linearize the curve y = f(x) by
                    by chosing b1 = b2 and n = 0. Here  f is one dimensional 
                    BB-form. 
----------------------------------------------------------------------------*/
void
Display2DFunction(a1,a2,b1,b2,m,n,weight,scaler,detx,dety,polygons)
double   weight[];               /* the coefficients of Bernstein polynomial*/
double a1,a2,b1,b2;              /* a1 <= x <= a2; b1 <= y <= b2;           */
int    m,n;                      /* degree in x and y direction m,n <= 5    */
double scaler[];                 /* the mapping parameter                   */
double detx,dety;                /* x and y dirction sizes in display the
                                    surface z = f(x,y)                      */
Polygons *polygons;              /* polygons chain  created                 */
{
double dxweight[30], dyweight[30],detx1,dety1; 
int   i,j,k,mx,my;
double   x,y,z,d,aa1,aa2,bb1,bb2;
Vertex  *vertices; 
Normal  *normals;
Polygons end,temp;

if (*polygons != NULL) {
   fprintf(stderr,
        "polygons need to be NUll pointer, now it is force to be NULL!\n");
   *polygons = NULL;
}
mx = (a2 - a1)/detx + 0.5;
my = (b2 - b1)/dety + 0.5;
if (mx == 0) mx = 1;
if (my == 0) my = 1;

detx1 = (a2 - a1)/mx;
dety1 = (b2 - b1)/my;

PartialxyBernstein2D(weight,m,n,a1,a2,b1,b2,dxweight,0);
PartialxyBernstein2D(weight,m,n,a1,a2,b1,b2,dyweight,1);


for (k = 0; k < mx; k++)         /* k loop     */
for (j = 0; j < my; j++){        /* j loop     */
   aa1 = a1 + k*detx1;  aa2 = aa1 + detx1;
   bb1 = b1 + j*dety1;  bb2 = bb1 + dety1;
   vertices = (Vertex *) malloc (4 * sizeof (Vertex));
   normals = (Normal *) malloc (4 * sizeof (Normal));

   x = (aa1-a1)/(a2 - a1);   y = (bb1 - b1)/(b2 - b1);
   vertices[0][0] = scaler[0]*aa1 + scaler[1];
   vertices[0][1] = scaler[0]*bb1 + scaler[2];
   vertices[0][2] = scaler[0]*BernsteinValue2D(weight,x,y,m,n) + scaler[3];
   normals[0][0] = BernsteinValue2D(dxweight,x,y,m-1,n);
   normals[0][1] = BernsteinValue2D(dyweight,x,y,m,n-1);
   normals[0][2] = -1.0;

    x = (aa1-a1)/(a2 - a1);   y = (bb2 - b1)/(b2 - b1);
   vertices[1][0] = scaler[0]*aa1 + scaler[1];
   vertices[1][1] = scaler[0]*bb2 + scaler[2];
   vertices[1][2] = scaler[0]*BernsteinValue2D(weight,x,y,m,n) + scaler[3];
   normals[1][0] = BernsteinValue2D(dxweight,x,y,m-1,n);
   normals[1][1] = BernsteinValue2D(dyweight,x,y,m,n-1);
   normals[1][2] = -1.0;

    x = (aa2-a1)/(a2 - a1);   y = (bb2 - b1)/(b2 - b1);
   vertices[2][0] = scaler[0]*aa2 + scaler[1];
   vertices[2][1] = scaler[0]*bb2 + scaler[2];
   vertices[2][2] = scaler[0]*BernsteinValue2D(weight,x,y,m,n) + scaler[3];
   normals[2][0] = BernsteinValue2D(dxweight,x,y,m-1,n);
   normals[2][1] = BernsteinValue2D(dyweight,x,y,m,n-1);
   normals[2][2] = -1.0;

   x = (aa2-a1)/(a2 - a1);   y = (bb1 - b1)/(b2 - b1);
   vertices[3][0] = scaler[0]*aa2 + scaler[1];
   vertices[3][1] = scaler[0]*bb1 + scaler[2];
   vertices[3][2] = scaler[0]*BernsteinValue2D(weight,x,y,m,n) + scaler[3];
   normals[3][0] = BernsteinValue2D(dxweight,x,y,m-1,n);
   normals[3][1] = BernsteinValue2D(dyweight,x,y,m,n-1);
   normals[3][2] = -1.0;
   x = 0.0;  y = 0.0;  z = 0.0;  d = 0.0;
   for (i = 0; i < 3 ; i++){
       x = x + normals[0][i]*normals[0][i];
       y = y + normals[1][i]*normals[1][i];
       z = z + normals[2][i]*normals[2][i];
       d = d + normals[3][i]*normals[3][i];
   }
   x = sqrt(x);  y = sqrt(y);   z = sqrt(z);   d = sqrt(d);
   for (i = 0; i < 3 ; i++){
       normals[0][i] = normals[0][i]/x;
       normals[1][i] = normals[1][i]/y;
       normals[2][i] = normals[2][i]/z;
       normals[3][i] = normals[3][i]/d;
   }
   if (b2 - b1 < 0.0002) {
      Vertex  *vertices1;
      vertices1 = (Vertex *) malloc (2 * sizeof (Vertex));
	  vertices[0][1] = 0.0;
	  vertices[3][1] = 0.0;
      for (i = 0; i < 3 ; i++){
          vertices1[0][i] = vertices[0][i];
          vertices1[1][i] = vertices[3][i];
      }
/*
      LineAbs(vertices[0][0],0.0,vertices[0][2],vertices[3][0],0.0,vertices[3][2]);
*/
      temp = CreatPolygons(2, vertices1, NULL);
      if (*polygons == NULL){
         *polygons = temp;
         end = temp;
      }  else {
         end->next = temp;
         end = temp;
      }
   }  else {
	  temp = CreatPolygons(4, vertices, normals);
      if (*polygons == NULL) {
		 *polygons = temp;
         end = temp;
      }  else {
	     end->next = temp;
         end = temp;
      }
	 /*
      DrawPolygon(4, vertices, normals);
	  */
   }
}       /* end j loop and k loop  */
}

/*----------------------------------------------------------------------------
CreatPolygons -- creat a polygon
----------------------------------------------------------------------------*/
Polygons 
CreatPolygons(n,vertices,normals)
int      n;                      /* the number of vertices, n > 0           */
Vertex   *vertices;              /* the vertices of the polygon being 
					                created                                 */
Normal   *normals;               /* the normals of the polygon if n = 1 or 2,
									the normal can be NULL                  */
{
Polygons result;

result = (Polygons)malloc(sizeof(struct _Polygons_));
result->numpts = n;
result->vertices = vertices;
result->normals = normals;
result->next = NULL;

return(result);
}

/*----------------------------------------------------------------------------
FreePolygons -- free  a polygon chain
----------------------------------------------------------------------------*/
void
FreePolygons(polygons)
Polygons polygons;               /* the given polygons                      */
{

Polygons next_polygons;
while (polygons != NULL) {
    if (polygons->vertices != NULL)   free(polygons->vertices);
    if (polygons->normals  != NULL)   free(polygons->normals);
	next_polygons = polygons->next;
	free(polygons);
	polygons = next_polygons;
}
}
