#ifndef MESH_H
#define MESH_H

#include <sys/types.h>
#include <Decimation/Decim/ct.h>
#ifdef __cplusplus
extern "C" {
#endif

#define VT_MAGIC 0x5ca4ef2b

/* supported mesh types */
typedef enum {
  /* 2d points & lines */
  VT_MESH_SCAT2D = 3,
  VT_MESH_LINE2D = 8,

  /* 2d meshes */
  VT_MESH_RECT2D = 0,
  VT_MESH_CURV2D = 1,
  VT_MESH_UNST2D = 2,
  VT_MESH_QUNST2D = 12,

  /* 3d points & lines */
  VT_MESH_SCAT3D = 7,
  VT_MESH_LINE3D = 9,

  /* surfaces in 3D */
  VT_MESH_RECT3D2D = 10,
  VT_MESH_SURF3D = 11,

  /* volume 3d meshes */
  VT_MESH_RECT3D = 4,
  VT_MESH_CURV3D = 5,
  VT_MESH_UNST3D = 6,

  /* remote meshes (through client/server connection) */
  VT_MESH_REMOTE = 99
} VTmeshType;

/* supported mesh types */
typedef enum {
  VT_UCHAR = 0,
  VT_CHAR = 1,
  VT_SHORT = 2,
  VT_LONG = 3,
  VT_FLOAT = 4,
  VT_COLOR = 5
} VTvarType;

typedef struct {
  union {
    u_char *ucdata;
    char *cdata;
    short *sdata;
    long *ldata;
    float *fdata;
    u_char (*rgbdata)[3];
  } data;
} VTline2Ddata, *VTline2DdataP;

typedef struct {
  double x, y;
} VTline2Dvdata, *VTline2DvdataP;

typedef struct {
  int vert[2];
} VTline2Dldata, *VTline2DldataP;

typedef struct {
  int nverts, nlines;
  VTline2DvdataP verts;
  VTline2DldataP lines;
  VTline2DdataP *data;
} VTline2Dtdata, *VTline2DtdataP;

/* 2D line data */
typedef struct {
  VTline2DtdataP *data; /* data for a particular timestep */
} VTline2Dmesh, *VTline2DmeshP;

/* data for a single variable, single timestep */
typedef struct {
  union {
    u_char *ucdata;
    char *cdata;
    short *sdata;
    long *ldata;
    float *fdata;
    u_char (*rgbdata)[3];
  } data;
} VTrect2Ddata, *VTrect2DdataP;

typedef struct {
  VTrect2DdataP *data; /* data for each variable, one timestep */
} VTrect2Dtdata, *VTrect2DtdataP;

/* 2D rectilinear data */
typedef struct {
  int xdim, ydim;              /* x,y dimensions */
  double xo, yo, xskip, yskip; /* x,y origin (bottom left) and strides */
  VTrect2DtdataP *data;        /* data for a particular timestep */
} VTrect2Dmesh, *VTrect2DmeshP;

/* 2D curvilinear data */
typedef struct {
  int xdim, ydim;
} VTcurv2Dmesh, *VTcurv2DmeshP;

typedef struct {
  double x, y;
} VTunst2Dvdata, *VTunst2DvdataP;

typedef struct {
  int vert[3];
  int adjtri[3];
} VTunst2Dtdata, *VTunst2DtdataP;

typedef struct {
  union {
    u_char *ucdata;
    char *cdata;
    short *sdata;
    long *ldata;
    float *fdata;
    u_char (*rgbdata)[3];
  } data;
} VTunst2Ddata, *VTunst2DdataP;

typedef struct {
  int nvert, ntri;
  /* mesh info */
  VTunst2DvdataP verts;
  VTunst2DtdataP tri;
  VTunst2DdataP *data;
} VTunst2Dmdata, *VTunst2DmdataP;

/* 2D unstructured data */
typedef struct {
  VTunst2DmdataP *data;
} VTunst2Dmesh, *VTunst2DmeshP;

typedef struct {
  double x, y;
} VTqunst2Dvdata, *VTqunst2DvdataP;

typedef struct {
  int vert[3];
  int adj[3];
} VTqunst2Dtdata, *VTqunst2DtdataP;

typedef struct {
  int vert[4];
  int adj[4];
} VTqunst2Dqdata, *VTqunst2DqdataP;

typedef struct {
  union {
    u_char *ucdata;
    char *cdata;
    short *sdata;
    long *ldata;
    float *fdata;
    u_char (*rgbdata)[3];
  } data;
} VTqunst2Ddata, *VTqunst2DdataP;

typedef struct {
  int nvert, ntri, nquad;
  /* mesh info */
  VTqunst2DvdataP verts;
  VTqunst2DtdataP tri;
  VTqunst2DqdataP quad;
  VTqunst2DdataP *data;
} VTqunst2Dmdata, *VTqunst2DmdataP;

/* 2D quad/tri unstructured data */
typedef struct {
  VTqunst2DmdataP *data;
} VTqunst2Dmesh, *VTqunst2DmeshP;

typedef struct {
  double x, y;
} VTscat2Dvdata, *VTscat2DvdataP;

typedef struct {
  union {
    u_char *ucdata;
    char *cdata;
    short *sdata;
    long *ldata;
    float *fdata;
    u_char (*rgbdata)[3];
  } data;
} VTscat2Ddata, *VTscat2DdataP;

typedef struct {
  VTscat2DdataP *data; /* data for n variables */
} VTscat2Dmdata, *VTscat2DmdataP;

/* 2D scattered data */
typedef struct {
  int nvert;
  VTscat2DvdataP verts;
  VTscat2DmdataP *data; /* data for n timesteps */
} VTscat2Dmesh, *VTscat2DmeshP;

typedef struct {
  union {
    u_char *ucdata;
    char *cdata;
    short *sdata;
    long *ldata;
    float *fdata;
    u_char (*rgbdata)[3];
  } data;
} VTline3Ddata, *VTline3DdataP;

typedef struct {
  double pt[3];
} VTline3Dvdata, *VTline3DvdataP;

typedef struct {
  int vert[2];
} VTline3Dldata, *VTline3DldataP;

typedef struct {
  int nverts, nlines;
  VTline3DvdataP verts;
  VTline3DldataP lines;
  VTline3DdataP *data;
} VTline3Dtdata, *VTline3DtdataP;

/* 3D line data */
typedef struct {
  VTline3DtdataP *data; /* data for a particular timestep */
} VTline3Dmesh, *VTline3DmeshP;

/* data for a single variable, single timestep */
typedef struct {
  union {
    u_char *ucdata;
    char *cdata;
    short *sdata;
    long *ldata;
    float *fdata;
    u_char (*rgbdata)[3];
  } data;
} VTrect3D2Ddata, *VTrect3D2DdataP;

typedef struct {
  VTrect3D2DdataP *data; /* data for each variable, one timestep */
} VTrect3D2Dtdata, *VTrect3D2DtdataP;

/* 3D2D rectilinear data */
typedef struct {
  int xdim, ydim; /* x,y dimensions */
  double xo, yo, zo, xskip[3], yskip[3];
  /* x,y origin (bottom left) and strides */
  VTrect3D2DtdataP *data; /* data for a particular timestep */
} VTrect3D2Dmesh, *VTrect3D2DmeshP;

/* data for a single variable, single timestep */
typedef struct {
  union {
    u_char *ucdata;
    char *cdata;
    short *sdata;
    long *ldata;
    float *fdata;
    u_char (*rgbdata)[3];
  } data;
} VTrect3Ddata, *VTrect3DdataP;

typedef struct {
  VTrect3DdataP *data; /* data for each variable, one timestep */
} VTrect3Dtdata, *VTrect3DtdataP;

/* 3D rectilinear data */
typedef struct {
  int xdim, ydim, zdim;
  double xo, yo, zo, xskip, yskip, zskip; /* x,y,z origin (bottom left) and
                                strides */
  VTrect3DtdataP *data;                   /* data for a particular timestep */
} VTrect3Dmesh, *VTrect3DmeshP;

/* 3D curvilinear data */
typedef struct {
  int xdim, ydim, zdim;
} VTcurv3Dmesh, *VTcurv3DmeshP;

typedef struct {
  double pt[3];
} VTsurf3Dvdata, *VTsurf3DvdataP;

typedef struct {
  int vert[3];
  int adjtri[3];
} VTsurf3Dtdata, *VTsurf3DtdataP;

typedef struct {
  union {
    u_char *ucdata;
    char *cdata;
    short *sdata;
    long *ldata;
    float *fdata;
    u_char (*rgbdata)[3];
  } data;
} VTsurf3Ddata, *VTsurf3DdataP;

typedef struct {
  int nvert, ntri;
  /* mesh info */
  VTsurf3DvdataP verts;
  VTsurf3DtdataP tri;
  VTsurf3DdataP *data;
} VTsurf3Dmdata, *VTsurf3DmdataP;

/* 3D surface data */
typedef struct {
  VTsurf3DmdataP *data;
} VTsurf3Dmesh, *VTsurf3DmeshP;

typedef struct {
  double pt[3];
} VTunst3Dvdata, *VTunst3DvdataP;

typedef struct {
  int vert[4];
  int adj[4];
} VTunst3Dtdata, *VTunst3DtdataP;

typedef struct {
  union {
    u_char *ucdata;
    char *cdata;
    short *sdata;
    long *ldata;
    float *fdata;
    u_char (*rgbdata)[3];
  } data;
} VTunst3Ddata, *VTunst3DdataP;

/* 3D unstructured data */
typedef struct {
  int nvert, ntet;
  /* mesh info */
  VTunst3DvdataP verts;
  VTunst3DtdataP tet;
  VTunst3DdataP *data;
} VTunst3Dmdata, *VTunst3DmdataP;

/* 3D unstructured data */
typedef struct {
  VTunst3DmdataP *data;
} VTunst3Dmesh, *VTunst3DmeshP;

/* 3D scattered data */
typedef struct {
  int npts;
} VTscat3Dmesh, *VTscat3DmeshP;

/* remote mesh */
typedef struct {
  int client;
  int meshnum;
} VTremotemesh, *VTremotemeshP;

typedef struct {
  char name[10]; /* name of problem variable */
  double min, max;
  VTvarType vtype; /* type of this variable */
} VTvar, *VTvarP;

/* a mesh */
typedef struct {
  VTmeshType mesh_type;
  double wc[3];     /* world center of mesh */
  double extent;    /* extent of mesh */
  int nvars, ntime; /* # of variables and timesteps */
  VTvarP vars;      /* variable info */
  double *times;    /* times */
  union {
    VTline2DmeshP line2d;
    VTscat2DmeshP scat2d;

    VTrect2DmeshP rect2d;
    VTcurv2DmeshP curv2d;
    VTunst2DmeshP unst2d;
    VTqunst2DmeshP qunst2d;

    VTline3DmeshP line3d;
    VTscat3DmeshP scat3d;

    VTrect3D2DmeshP rect3d2d;
    VTsurf3DmeshP surf3d;

    VTrect3DmeshP rect3d;
    VTcurv3DmeshP curv3d;
    VTunst3DmeshP unst3d;

    VTremotemeshP remote;
  } mesh;
} VTmesh, *VTmeshP;

int VTwriteMesh(char *, VTmeshP);
int VTwriteMeshTime(char *, VTmeshP, int);
VTmeshP VTreadMesh(char *);
VTmeshP VTreadPBM(char *);
VTmeshP VTreadDEM(char *);
VTmeshP VTreadSLC(char *);
VTmeshP VTslice2mesh(CTSlice);
VTmeshP VTvol2mesh(CTVolume);
int VTgenTStrips(VTmeshP, int);
void VTmeshMinMaxVars(VTmeshP);

/* access to mesh independent data */
#define VTmeshMeshType(vtmp) ((vtmp)->mesh_type)
#define VTmeshExtent(vtmp) ((vtmp)->extent)
#define VTmeshNvars(vtmp) ((vtmp)->nvars)
#define VTmeshVarType(vtmp, v) ((vtmp)->vars[v].vtype)
#define VTmeshVarMin(vtmp, v) ((vtmp)->vars[v].min)
#define VTmeshVarMax(vtmp, v) ((vtmp)->vars[v].max)
#define VTmeshVarName(vtmp, v) ((vtmp)->vars[v].name)

#define VTmeshNtime(vtmp) ((vtmp)->ntime)
#define VTmeshTime(vtmp, t) ((vtmp)->times[t])

/* access to mesh specific data */
/*
#define VTmeshType(vtmp)      ((vtmp)->mesh_type)
*/
#define VTmeshLine2D(vtmp) ((vtmp)->mesh.line2d)
#define VTmeshRect2D(vtmp) ((vtmp)->mesh.rect2d)
#define VTmeshCurv2D(vtmp) ((vtmp)->mesh.curv2d)
#define VTmeshUnst2D(vtmp) ((vtmp)->mesh.unst2d)
#define VTmeshScat2D(vtmp) ((vtmp)->mesh.scat2d)
#define VTmeshLine3D(vtmp) ((vtmp)->mesh.line3d)
#define VTmeshRect3D(vtmp) ((vtmp)->mesh.rect3d)
#define VTmeshCurv3D(vtmp) ((vtmp)->mesh.curv3d)
#define VTmeshUnst3D(vtmp) ((vtmp)->mesh.unst3d)
#define VTmeshScat3D(vtmp) ((vtmp)->mesh.scat3d)

#ifdef __cplusplus
};
#endif

#endif
