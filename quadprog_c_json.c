/*
  gcc -o quadprog quadprog.c && ./quadprog *.json
 */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "jsmn.h"
#include "qpgen2_.h"

static int verbose = 0;
static int gidx = 0;
static char *gfmt[] = {"%s%g", "%s%15.13g", "%s%20.18g"};

double epsilon = 0.0;

typedef struct {

  // problem
  int n; // variables
  int m; // constraints
  int meq; // equal constraints
  double *G; // n x n
  double *a; // n x 1
  double *C; // n x m
  double *b; // m x 1
  int factorized; // falg

  // solution
  double *opt; // n x 1, optimal variables
  double value; // optimal value
  double *unc; // n x 1, unconstrained solution
  double *l; // n x 1, lagrangian
  int *iter; // 2, iteration data
  int niact; // #active constraints
  int *iact; // < n x 1, active constraint indexes
} qp;

typedef struct {
  char *chars;
  int nchar;
  jsmntok_t *tokens;
  int ntoken;
} json;

json *json_free(json *js) {
  if (js) {
    if (js->chars) free(js->chars);
    if (js->tokens) free(js->tokens);
    free(js);
  }
  return NULL;
}

json *json_from_filename(char *fn) {
  FILE *f;
  jsmn_parser p;
  jsmntok_t *t;
  int r;
  int ln;
  char *str;
  json *js;

  js = malloc(sizeof(*js));
  if (!js)
    return NULL;

  js->ntoken = JSMN_ERROR_NOMEM;
  js->chars = NULL;
  js->nchar = 0;
  js->tokens = NULL;
  js->ntoken = 0;
  
  f = fopen(fn, "rb");
  if (f == NULL)
    return js;
  fseek(f, 0, SEEK_END);
  ln = ftell(f);
  fseek(f, 0, SEEK_SET);
  if ( (str = malloc(ln+1)) == NULL) {
    fclose(f);
    return js;
  }

  js->chars = str;
  js->nchar = ln;
  if (fread(str, 1, ln, f) != ln) {
    fclose(f);
    return js;
  }
  str[ln] = '\0';
  fclose(f);
  
  jsmn_init(&p);
  r = jsmn_parse(&p, str, ln, NULL, 0);
  switch(r) {
  case JSMN_ERROR_INVAL:
  case JSMN_ERROR_NOMEM:
  case JSMN_ERROR_PART :
    return js;
  default: break;
  }

  if ( (t = (jsmntok_t *)malloc(sizeof(*t) * r)) == NULL)
    return js;
  js->tokens = t;
  
  jsmn_init(&p);
  r = jsmn_parse(&p, str, ln, t, r);
  switch (r) {
  case JSMN_ERROR_INVAL:
  case JSMN_ERROR_NOMEM:
  case JSMN_ERROR_PART :
    return js;
  default: break;
  }
  js->ntoken = r;

  return js;
}

void token_info(char *l, char *p, jsmntok_t *t) {
  int sz = t->end - t->start;
  p += t->start;

  printf("INFO: %s ", l);
  while (sz--)
    printf("%c", *p++);
  printf("\n");
}

jsmntok_t *findinobject(char *name, json *js) {
  jsmntok_t *t = js->tokens;
  jsmntok_t *lt = &t[js->ntoken-1];
  long nln = strlen(name);
  char *str = js->chars;
  long ln;
  int et;

  if (t++->type != JSMN_OBJECT || js->ntoken <= 0)
    return NULL;

  while (t <= lt) {
    if (t->type != JSMN_STRING)
      return NULL;

    ln = t->end - t->start;
    if (ln == nln && strncmp(str+t->start, name, ln) == 0)
      return t+1;

    t++;
    et = t->end;
    t++;

    while (t->end < et)
      t++;
  }

  return NULL;
}

void getdarray(jsmntok_t *t, json *js, double *v) {
  int i;
  char *str;

  if (!t)
    return;

  if (t->type != JSMN_ARRAY)
    return;

  str = js->chars;

  for (i = 1; i <= t->size; i++) {
    if (t[i].type != JSMN_PRIMITIVE)
      return;

    switch(str[t[i].start]) {
    case 't': case 'T': v[i-1] = 1; break;
    case 'f': case 'F': v[i-1] = 0; break;
    case '-': case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9': case '.': case '+':
      v[i-1] = strtod(str+t[i].start, NULL);
      break;
    }
  }
}

void getiarray(jsmntok_t *t, json *js, int *v) {
  int i;
  char *str;

  if (!t)
    return;

  if (t->type != JSMN_ARRAY)
    return;

  str = js->chars;

  for (i = 1; i <= t->size; i++) {
    if (t[i].type != JSMN_PRIMITIVE)
      return;

    switch(str[t[i].start]) {
    case 't': case 'T': v[i-1] = 1; break;
    case 'f': case 'F': v[i-1] = 0; break;
    case '-': case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9': case '.': case '+':
      v[i-1] = strtol(str+t[i].start, NULL, 0);
      break;
    }
  }
}

int geti0(char *name, json *js, int def) {
  int i = def;
  char *str = js->chars;

  jsmntok_t *t = findinobject(name, js);
  if (t[0].type == JSMN_PRIMITIVE) {
    switch(str[t[0].start]) {
    case 't': case 'T': i = 1; break;
    case 'f': case 'F': i = 0; break;
    case '-': case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9': case '.': case '+':
      i = strtol(str+t[0].start, NULL, 0);
      break;
    }
  }

  return i;
}

double getd0(char *name, json *js, double def) {
  double d = def;
  char *str = js->chars;

  jsmntok_t *t = findinobject(name, js);
  if (t[0].type == JSMN_PRIMITIVE) {
    switch(str[t[0].start]) {
    case 't': case 'T': d = 1.0; break;
    case 'f': case 'F': d = 0.0; break;
    case '-': case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9': case '.': case '+':
      d = strtod(str+t[0].start, NULL);
      break;
    }
  }

  return d;
}

double *getd1(char *name, json *js, int *sz) {
  jsmntok_t *t = findinobject(name, js);
  double *v = NULL;
  
  *sz = 0;
  if (!t || t->size == 0)
    return v;
  v = malloc(sizeof(*v)*t->size);
  if (!v)
    return v;

  getdarray(t, js, v);
  *sz = t->size;
  return v;
}

int *geti1(char *name, json *js, int *sz) {
  jsmntok_t *t = findinobject(name, js);
  int *v = NULL;
  
  *sz = 0;
  if (!t || t->size == 0)
    return v;
  v = malloc(sizeof(*v)*t->size);
  if (!v)
    return v;

  getiarray(t, js, v);
  *sz = t->size;
  return v;
}

double *getd2(char *name, json *js, int *sz1, int *sz2) {
  double *v = NULL;
  double *w;
  double *a;
  jsmntok_t *t = findinobject(name, js);
  char *str = js->chars;
  jsmntok_t *lt = &t[js->ntoken-1];
  int n1, n2;
  int i;
  int j;

  *sz1 = 0;
  *sz2 = 0;

  if (!t || t+3 > lt)
    return v;

  if (t[0].type != JSMN_ARRAY ||
      t[1].type != JSMN_ARRAY)
    return v;

  n1 = t[0].size;
  n2 = t[1].size;
  if ( (w = malloc(sizeof(*v) * n1 * n2)) == 0 )
    return v;
  if ( (a = malloc(sizeof(*a) * n2)) == 0 )
    return v;
  v = w;

  t++;
  for (i = 0; i < n1; i++) {
    if (t->size != n2)
      return v;
    getdarray(t, js, a);
    for (j = 0; j < n2; j++)
      v[i + j * n1] = a[j];
    t += n2 + 1;
  }
  *sz1 = n1;
  *sz2 = n2;

  free(a);
  return v;
}

void putd0(char *l, double v) {
  int i;

  printf("INFO: ");
  printf(gfmt[gidx], l, v);
  printf("\n");
}

void putd1(char *l, double *v, int n) {
  int i;

  if (!v)
    return;

  printf("INFO: %s: [", l);
  for (i = 0; i < n; i++)
    printf(gfmt[gidx], i ? ", " : "", v[i]);
  printf("]\n");
}

void putb0(char *l, int v) {
  int i;

  printf("INFO: %s: %s\n", l, v ? "true" : "false");
}

void puti0(char *l, int v) {
  int i;

  printf("INFO: %s: %d\n", l, v);
}

void puti1(char *l, int *v, int n) {
  int i;

  if (!v)
    return;

  printf("INFO: %s: [", l);
  for (i = 0; i < n; i++)
    printf("%s%d", (i?", ":""),v[i]);
  printf("]\n");
}

void putd2(char *l, double *v, int n1, int n2) {
  int i;
  int j;

  if (!v)
    return;

  printf("INFO: %s: [", l);
  for (i = 0; i < n1; i++) {
    printf("%s[", i ? ", ": "");
    for (j = 0; j < n2; j++)
      printf(gfmt[gidx], j ? ", " : "", v[ i + j * n1]);
    printf("]");
  }
  printf("]\n");
}

  
qp *qp_free(qp *p) {
  if (p) {
    if (p->G) free(p->G);
    if (p->a) free(p->a);
    if (p->C) free(p->C);
    if (p->b) free(p->b);
    if (p->opt) free(p->opt);
    if (p->unc) free(p->unc);
    if (p->l) free(p->l);
    if (p->iter) free(p->iter);
    if (p->iact) free(p->iact);
    free(p);
  }
  return NULL;
}

qp *qp_from_json(json *js) {
  qp *p = calloc(1, sizeof(*p));
  int v;
  if (!p)
    return p;
  p->G = getd2("G", js, &p->n, &v);
  if (p->n != v) { // G must be square
    printf("FAIL: G not square\n");
    return qp_free(p);
  }
  p->a = getd1("a", js, &v);
  if (p->n != v) { // a must align with G
    printf("FAIL: a not like G\n");
    return qp_free(p);
  }
  p->C = getd2("C", js, &v, &p->m);
  if (p->C) {
    if (p->n != v) { // C must align with G
      printf("FAIL: C not like G\n");
      return qp_free(p);
    }
  } else
    p->m = 0;
  p->b = getd1("b", js, &v);
  if (p->b && p->m != v) { // b must align with C
    printf("FAIL: b not like C\n");
    return qp_free(p);
  }
  p->meq = geti0("meq", js, 0);
  p->factorized = geti0("factorized", js, 0);

  p->opt = getd1("solution", js, &v);
  if (p->opt && p->n != v) { //optional
    printf("FAIL: opt not like G\n");
    return qp_free(p);
  }

  p->value = getd0("value", js, NAN);
  p->unc = getd1("unconstrained.solution", js, &v);
  if (p->unc && p->n != v) { //optional
    printf("FAIL: unc not like G\n");
    return qp_free(p);
  }
  p->l = getd1("Lagrangian", js, &v);
  if (p->l && p->m != v) { //optional
    printf("FAIL: l not like b\n");
    return qp_free(p);
  }
  p->iter = geti1("iterations", js, &v);
  if (p->iter && v != 2) { //optional
    printf("FAIL: iter not 2\n");
    return qp_free(p);
  }
  p->iact = geti1("iact", js, &p->niact);
  if (p->iact && p->niact > p->m) { //optional
    printf("FAIL: iact too big\n");
    return qp_free(p);
  }

  return p;
}

void qp_info(char *l, qp *p) {
  if (l)
    printf("INFO:\nINFO: %s:\n", l);

  if (p) {
    if (p->G) putd2("G", p->G, p->n, p->n);
    if (p->a) putd1("a", p->a, p->n);
    if (p->C) putd2("C", p->C, p->n, p->m);
    if (p->b) putd1("b", p->b, p->m);
    putd0("meq", p->meq);
    putb0("factorized", p->factorized);
    if (p->opt) putd1("opt", p->opt, p->n);
    putd0("value", p->value);
    if (p->unc) putd1("unc", p->unc, p->n);
    if (p->iter) puti1("iterations", p->iter, 2);
    if (p->l) putd1("l", p->l, p->m);
    if (p->iact) puti1("iact", p->iact, p->niact);
  }
}

void *memcopy(void *p, int ln) {
  void *q = malloc(ln);
  if (q)
    memcpy(q, p, ln);
  return q;
}

qp *qp_copy_problem(qp *p) {
  qp *q = calloc(1, sizeof(*q));
  if (q) {
    q->n = p->n;
    q->m = p->m;
    q->meq = p->meq;
    q->factorized = p->factorized;
    q->G = memcopy(p->G, p->n * p->n * sizeof(*p->G));
    q->a = memcopy(p->a, p->n * sizeof(*p->a));
    q->C = memcopy(p->C, p->n * p->m * sizeof(*p->C));
    q->b = memcopy(p->b, p->m * sizeof(*p->b));

    q->opt = calloc(q->n, sizeof(*q->opt));
    q->value = 0;
    q->unc = calloc(q->n, sizeof(*q->unc));
    q->l = calloc(p->m, sizeof(*q->l));
    q->iter = calloc(2, sizeof(*q->iter));
    q->niact = 0;
    q->iact = calloc(q->m, sizeof(*q->iact));
  }
  return q;
}

int samei1(int *p, int *q, int n) {
  int i;
  for (i = 0; i < n; i++, p++, q++) {
    if (*p != *q)
      return 0;
  }
  return 1;
}

int samed1(double *p, double *q, int n) {
  int i;
  double diff;
  double threshold;
  for (i = 0; i < n; i++, p++, q++) {
    diff = *p-*q;
    if (diff < 0) diff = -diff;
    threshold = epsilon + 1e-10 * (*q < 0 ? -*p : *p);
    if (diff > threshold)
      return 0;
  }
  return 1;
}

int qp_same(qp *p, qp *q, int *minor) {

  *minor = 0;

  if (verbose > 0) {
    putd1("p opt", p->opt, p->n); printf("\n");
    putd1("q opt", q->opt, q->n); printf("\n");
  }
  if (samed1(p->opt, q->opt, p->n) == 0)
    return 0;

  if (verbose > 0) {
    putd1("p unc", p->unc, p->n); printf("\n");
    putd1("q unc", q->unc, q->n); printf("\n");
  }
  if (samed1(p->unc, q->unc, p->n) == 0)
    return 0;

  if (verbose > 0) {
    putd0("p val", p->value);
    putd0("q val", q->value);
  }
  if (samed1(&p->value, &q->value, 1) == 0)
    return 0;

  if (verbose > 0) {
    putd1("p lag", p->l, p->m); printf("\n");
    putd1("q lag", q->l, q->m); printf("\n");
  }
  if (samed1(p->l, q->l, p->m) == 0)
    (*minor)++;
  

  if (verbose > 0) {
    puti1("p itr", p->iter, 2); printf("\n");
    puti1("q itr", q->iter, 2); printf("\n");
  }
  if (samei1(p->iter, q->iter, 2) == 0)
    (*minor)++;
  
  if (verbose > 0) {
    puti1("p act", p->iact, p->niact); printf("\n");
    puti1("q act", q->iact, q->niact); printf("\n");
  }
  if (samei1(p->iact, q->iact, p->niact) == 0)
    (*minor)++;

  return 1;
}

int qptest(char *fn, json *js, int *minor) {
  qp *p = qp_from_json(js);

  if (!p) {
    *minor = -1;
    return 0;
  }

  qp *q;
  int tst = 1;
  int r = p->n > p->m ? p->m : p->n;
  double fake;
  double *work = malloc((2*p->n + 2*p->m + r*(r+5)/2) * sizeof(*work));
  int i;

  if (verbose > 1)
    qp_info("problem", p);
  
  q = qp_copy_problem(p);

  if (!q->C) q->C = &fake;
  if (!q->b) q->b = &fake;
  
  for (i = 0; i < q->n; i++)
    q->unc[i] = q->a[i];

  qpgen2_(q->G, q->unc, q->n, q->opt, q->l, &q->value,
	  q->C, q->b, q->m, q->meq,
	  q->iact, &q->niact, q->iter,
	  work, p->factorized);

  if (q->C == &fake) q->C = NULL;
  if (q->b == &fake) q->b = NULL;

  if (verbose > 1)
    qp_info("solution", q);

  free(work);

  tst = qp_same(p, q, minor);
  qp_free(p);
  qp_free(q);
  return tst;
}

void setoption(char *prg, char *opt) {
  switch(opt[1]) {
  case 'f':
    switch(opt[2]) {
    case 0: gidx = 1; break;
    case '0': gidx = 0; break;
    case '1': gidx = 1; break;
    case '2': gidx = 2; break;
    default : printf("INFO: ignoring option %s\n", opt); break;
    }
    break;

  case 'v':
    switch(opt[2]) {
    case 0: verbose = 1; break;
    case '0': verbose = 0; break;
    case '1': verbose = 1; break;
    case '2': verbose = 2; break;
    case '3': verbose = 3; break;
    default : printf("INFO: ignoring option %s\n", opt); break;
    }
    break;

  case 'h':
    printf("INFO: %s [-h] [-v#] [-f#] file.json ...\n", prg);
    printf("INFO: Finish\n");
    exit(1);

  default:
    printf("INFO: ignoring option %s\n", opt);
    break;
  }
}

int main(int argc, char *argv[]) {
  char *prg = argv[0];
  char *str;
  long ln;
  json *r;
  jsmntok_t *t;
  int minor;
  int options = 1;
  int fail = 0;

  printf("INFO: Start\n");

  epsilon = calculate_vsmall();
  while (argv++, --argc) {

    if (options && argv[0][0] == '-') {
      setoption(prg, argv[0]);
      continue;
    } else
      options = 0;

    if ( (r = json_from_filename(*argv)) == NULL) {
      printf("SKIP: %s - cannot read\n", *argv);
      continue;
    }

    if (r->ntoken <= 0) {
      printf("SKIP: %s - cannot parse\n", *argv);
      r = json_free(r);
      continue;
    }

    if (verbose) {
      if ( (t = findinobject("source", r)) != NULL)
	token_info(*argv, r->chars, t);

      if ( (t = findinobject("notes", r)) != NULL)
	token_info(*argv, r->chars, t);
    }

    minor = 0;
    if (qptest(*argv, r, &minor) == 0) {
      if (minor == -1)
	printf("SKIP: %s - bad configuration\n", *argv);
      else {
	printf("FAIL: %s\n", *argv);
	fail = 1;
      }
      break;
    }

    printf("PASS: %s", *argv);
    if (minor > 0)
      printf(", %d minor difference%s", minor, minor == 1 ? "" : "s");
    printf("\n");

    r = json_free(r);
  }

  printf("INFO: Finish: %s\n", fail ? "FAILED" : "SUCCESS");

  exit(fail);
}
