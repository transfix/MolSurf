/* defines for compatibility with 32-bit C-compilers */

#ifndef min
#define min(x, y) ((x) <= (y) ? (x) : (y))
#endif
#ifndef max
#define max(x, y) ((x) >= (y) ? (x) : (y))
#endif

#if !defined(__TURBOC__) && !defined(NDPC)
#define bioskey(x)
#endif

#ifndef __TURBOC__
#define far
#define near
#define huge
#define farmalloc malloc
#define farfree free
#define farcoreleft coreleft

#ifdef NDPC
#define cprintf __cprintf
#define coreleft() 0 /* not available */

/* colors */
enum COLORS {
  BLACK, /* dark colors */
  BLUE,
  GREEN,
  CYAN,
  RED,
  MAGENTA,
  BROWN,
  LIGHTGRAY,
  DARKGRAY, /* light colors */
  LIGHTBLUE,
  LIGHTGREEN,
  LIGHTCYAN,
  LIGHTRED,
  LIGHTMAGENTA,
  YELLOW,
  WHITE
};
extern int __k; /* scratch var. defined by NDP C */
#define setcolor set_color
#define lineto draw
#define moveto move
#define line(x0, y0, x1, y1)                                                   \
  {                                                                            \
    move(x0, y0);                                                              \
    draw(x1, y1);                                                              \
  }
#define clearviewport clear
#define closegraph text_mode

int stricmp(const char *s0, const char *s1);
int cprintf(const char *format, ...);
void gotoxy(int x, int y);
unsigned short bioskey(int mode);
#endif /* #ifdef NDPC */
#endif /* #ifndef __TURBOC__ */

/************************************************************************
 * EPSILON should be the smallest a such that 1.0+a!=1.0. (see epsilon.c)
 *************************************************************************/
#define EPSILON 1e-15
