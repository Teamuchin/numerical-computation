# Numerical Computation

Python implementations of two classical numerical methods: cubic spline interpolation and least-squares curve fitting.

**Course:** CENG216 — Numerical Computation  
**Institution:** İzmir Institute of Technology (IYTE) — İzmir, Türkiye

## Files

| File | Method |
|---|---|
| `natural-cubic-spline.py` | Builds a natural cubic spline through a set of points. Solves for the second derivatives by assembling the tridiagonal system and reducing it with `numpy.linalg.solve`, then derives the per-segment `a, b, c, d` coefficients and plots each segment |
| `func.py` | Shared helpers used by the spline construction |
| `piecewise-spline.py` | Evaluates and plots an explicit piecewise-cubic spline from its known segment polynomials, to check continuity at the knots |
| `least-squares-fit.py` | Least-squares straight-line fit. Forms the normal equations **XtX b = Xty** and solves them directly with `numpy.linalg.solve`, then plots the fit against the data |

Both spline scripts plot with Matplotlib so the interpolation can be inspected visually.

## Running

```bash
pip install numpy matplotlib
python natural-cubic-spline.py
```

The methods follow *Numerical Analysis* (Sauer), the course reference.

---

Submitted reports, worksheets and lecture material are archived outside this
repository rather than committed, so the repo stays code-only.
