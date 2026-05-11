/*  The purpose.  This lib contains the functions of linear approximations
                     of the iso-cuvre of a 2D BB-form polynomial and  iso-
                     surface of a 3D BB-form polynomial. The linear approx-
                     mation will be transformed by scalering and translatiton
                     in order to fit the screen.
                         Let P = (x,y,z) be a point in space. Then P is map
                     to
                         scaler[0]*x + scaler[1],
                         scaler[0]*y + scaler[2],
                         scaler[0]*z + scaler[3],
                     where scaler is a user chosen parameters, which will be
                     called mapping parameter. In most of the cases, one may
                     take
                            scaler[0] = 1.0; scaler[1] = 0.0;
                            scaler[2] = 0.0; scaler[3] = 0.0;
                     It means the transform is identical
*/

#ifndef DISPLAY_H
#define DISPLAY_H

typedef double Vertex[3]; /* the vertex of a polygon                */
typedef double Normal[3]; /* the normal at a vertex                 */

typedef struct _Polygons_ {
  int numpts;       /* the number of vertices                 */
  Vertex *vertices; /* the  vertices                          */
  Normal *normals;  /* the normals at the vertices            */
  struct _Polygons_ *next;
} *Polygons;

typedef struct PolygonLink {
  Polygons head;
  Polygons tail;
} *PolyLink;

void DisplayBernsL123(), /* Linearize a iso-surface of a BB-poly   */
    DisplayBerns123(),   /* Linearize a iso-surface of a BB-poly   */
    DisplayBerns2D123(), /* Linearize a iso-curve of a BB-poly     */
    Display2DFunction(), /* Linearize 2D function z=f(x,y), BB-poly*/
    FreePolygons(),      /* free a polygon chain                   */

    DisplayBerns123in(), DisplayBerns2D123in(), DisplayPolygons(),
    DisplayPolygons2D(),

    FindRoots2(), FindRoots2More(), FindRoots3(), FindRoots4(), BernsteinRoot(),
    BernsteinRootHigh(), FindRoots23Edges(), FindRoots2D123Edges(),
    Intersection12Edge(), Intersection12Edge1(), Intersection4Edge(),
    Intersection4Edge1(),

    Partition(), Partition2D(), GotPoint(),

    PartialBernstein(),     /* partial derivative in one direction    */
    PartialxyzBernstein(),  /* partial derivatives in three directions*/
    PartialxyBernstein2D(); /* partial derivatives in two directions  */

int PointsNum(), PointsNum2D(), which();

Polygons CreatPolygons(); /* creat a polygon                        */

#endif
