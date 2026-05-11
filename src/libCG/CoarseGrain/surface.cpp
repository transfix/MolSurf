#include <Utility/utility.h>
#include <libCG/CoarseGrain/geometry.h>
#include <libCG/CoarseGrain/surface.h>

using namespace GEOMETRY;

Surface::Surface() {}

Surface::~Surface() {
  for (int i = 0; i < m_Vertices.size(); i++) {
    delete m_Vertices[i];
  }
  for (int i = 0; i < m_Faces.size(); i++) {
    delete m_Faces[i];
  }
}

void Surface::addPoint(Vertex *point) {
  Vertex *p = new Vertex();
  p->Point[0] = point->Point[0];
  p->Point[1] = point->Point[1];
  p->Point[2] = point->Point[2];
  p->Normal[0] = point->Normal[0];
  p->Normal[1] = point->Normal[1];
  p->Normal[2] = point->Normal[2];
  m_Vertices.push_back(p);
}

void Surface::addFacet(Face *facet) {
  Face *f = new Face();
  f->Index[0] = facet->Index[0];
  f->Index[1] = facet->Index[1];
  f->Index[2] = facet->Index[2];
  f->AdjTri[0] = facet->AdjTri[0];
  f->AdjTri[1] = facet->AdjTri[1];
  f->AdjTri[2] = facet->AdjTri[2];
  f->whichnml[0] = facet->whichnml[0];
  f->whichnml[1] = facet->whichnml[1];
  f->whichnml[2] = facet->whichnml[2];
  m_Faces.push_back(f);
}
