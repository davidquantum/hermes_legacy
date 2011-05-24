#include "hermes2d.h"

using namespace WeakFormsH1;

/* Right-hand side */

class CustomExactFunction {
public:
  CustomExactFunction(double k, double alpha): k(k), alpha(alpha) { };

  double fn(double x, double y) {
    if (x <= 0) return cos(k * y);
    else return cos(k * y) + pow(x, alpha);
  }

  double k, alpha;
};

class CustomRightHandSide : public HermesFunction
{
public:
  CustomRightHandSide(double k, double alpha)
    : HermesFunction(), k(k), alpha(alpha) {
    cef = new CustomExactFunction(k, alpha);
  };

  virtual double value(double x, double y) const {
    if (x < 0) return cef->fn(x, y) * k * k;
    else return cef->fn(x, y) * k * k
                - alpha *(alpha - 1) * pow(x, alpha - 2.)
                - k * k * pow(x, alpha);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  ~CustomRightHandSide() { delete cef; }

  CustomExactFunction* cef;
  double k, alpha;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double k, double alpha)
    : ExactSolutionScalar(mesh), k(k), alpha(alpha) {
    cef = new CustomExactFunction(k, alpha);
  };

  virtual double value(double x, double y) const {
    return cef->fn(x, y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    if (x <= 0) dx = 0;
    else dx = alpha * pow(x, alpha - 1);
    dy = -sin(k * y) * k;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  ~CustomExactSolution() { delete cef; }

  CustomExactFunction* cef;
  double k, alpha;
};

