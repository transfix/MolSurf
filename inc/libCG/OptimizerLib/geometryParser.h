#ifndef GEOMETRYPARSER_H
#define GEOMETRYPARSER_H

#include <libCG/CoarseGrain/surface.h>
#include <stdio.h>

namespace GEOMETRY {
class Vertex;
class Face;
class Surface;
class GeometryParser {
public:
  GeometryParser();
  Surface *SurfaceTriangulation(FILE *fp);
  int ReadHowmanyComponents(FILE *fp);

private:
  int howmany, normalFlag;
};
} // namespace GEOMETRY

#endif // GEOMETRYPARSER_H
