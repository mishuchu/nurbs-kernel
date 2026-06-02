# 任务矩阵 — The NURBS Book 算法清单

## 章节汇总

| 章节 | 内容 | 算法数量 | 难度 |
|------|------|----------|------|
| Ch4 | B-spline 基数函数 | 5 | ⭐⭐ |
| Ch5 | NURBS 曲线 | 6 | ⭐⭐⭐ |
| Ch6 | NURBS 曲面 | 4 | ⭐⭐⭐⭐ |
| Ch7 | 曲线操作 | 4 | ⭐⭐⭐ |
| Ch8 | 曲面操作 | 3 | ⭐⭐⭐⭐ |
| Ch9 | 曲线曲面逼近 | 2 | ⭐⭐⭐ |
| Ch10 | 插值与逼近 | 3 | ⭐⭐⭐⭐ |
| Ch11 | 激光扫描数据处理 | 3 | ⭐⭐⭐ |
| Ch12 | 生产系统集成 | 3 | ⭐⭐⭐ |

**总计：33 个核心算法**

## 详细算法清单

### Ch4 — B-spline 基数函数
| 算法 | 函数签名 | 依赖 | 难度 |
|------|---------|------|------|
| B-样条基数函数 | `bspline_basis(u, p, U, k)` | 无 | ⭐ |
| 节点插入 | `insert_knot(u, U, p, Qw, k)` | Ch4 | ⭐⭐ |
| 节点细化 | `refine_knot_vector(U, p, X, a, b)` | Ch4 | ⭐⭐ |
| 节点删除 | `remove_knot(u, U, p, num, toler)` | Ch4 | ⭐⭐ |
| 升阶 | `elevate_degree(U, p, Qw, t)` | Ch4 | ⭐⭐ |

### Ch5 — NURBS 曲线
| 算法 | 函数签名 | 依赖 | 难度 |
|------|---------|------|------|
| 曲线构造 | `construct_curve(p, U, Pw)` | Ch4 | ⭐⭐ |
| 曲线求导 | `curve_derivs(u, p, U, Pw, d)` | Ch4 | ⭐⭐⭐ |
| 曲线节点插入 | `curve_knot_insertion(u, p, U, Qw)` | Ch4 | ⭐⭐ |
| 曲线升阶 | `curve_degree_elevation(p, U, Qw, t)` | Ch4 | ⭐⭐⭐ |
| 曲线反算 | `curve_inversion(u, p, U, Qw)` | Ch5构造 | ⭐⭐⭐ |
| 曲线积分 | `curve_integrate(p, U, Pw)` | Ch5求导 | ⭐⭐ |

### Ch6 — NURBS 曲面
| 算法 | 函数签名 | 依赖 | 难度 |
|------|---------|------|------|
| 曲面构造 | `construct_surface(m, n, U, V, PQw)` | Ch5 | ⭐⭐⭐ |
| 曲面求导 | `surface_derivs(u, v, p, q, U, V, PQw, d)` | Ch6 | ⭐⭐⭐⭐ |
| 曲面节点插入 | `surface_knot_insertion*(u_dir, v_dir, ...)` | Ch6 | ⭐⭐⭐ |
| 曲面升阶 | `surface_degree_elevation*(p, q, ...)` | Ch6 | ⭐⭐⭐⭐ |

### Ch7 — 曲线操作
| 算法 | 函数签名 | 依赖 | 难度 |
|------|---------|------|------|
| 曲线拆分 | `curve_subdivide(u, p, U, Qw)` | Ch5 | ⭐⭐ |
| 曲线合并 | `curve_merge(c1, c2)` | Ch7拆分 | ⭐⭐ |
| 重参数化 | `reparametrize(u_new, p, U, Qw)` | Ch5 | ⭐⭐⭐ |
| 曲线拟合 | `curve_fitting(points, p, k)` | Ch7 | ⭐⭐⭐ |

### Ch8 — 曲面操作
| 算法 | 函数签名 | 依赖 | 难度 |
|------|---------|------|------|
| 曲面分解 | `surface_decompose(p, q, U, V, PQw)` | Ch6 | ⭐⭐⭐⭐ |
| 曲面光顺 | `surface_smooth(p, q, U, V, PQw)` | Ch8分解 | ⭐⭐⭐⭐ |
| 曲面拼接 | `surface_join(s1, s2, dir)` | Ch8分解 | ⭐⭐⭐⭐ |

### Ch9 — 曲线曲面数据逼近
| 算法 | 函数签名 | 依赖 | 难度 |
|------|---------|------|------|
| 曲线逼近 | `curve_approximate(points, p, k)` | Ch7拟合 | ⭐⭐⭐ |
| 曲面逼近 | `surface_approximate(points, m, n)` | Ch9曲线 | ⭐⭐⭐⭐ |

### Ch10 — 插值与逼近
| 算法 | 函数签名 | 依赖 | 难度 |
|------|---------|------|------|
| 全局插值 | `global_interpolate(pts, p)` | Ch5 | ⭐⭐⭐ |
| Hermite插值 | `hermite_interpolate(P0, T0, P1, T1, p)` | Ch10全局 | ⭐⭐⭐ |
| 最小二乘逼近 | `least_squares_approx(points, p, n)` | Ch9 | ⭐⭐⭐⭐ |

### Ch11 — 激光扫描数据处理
| 算法 | 函数签名 | 依赖 | 难度 |
|------|---------|------|------|
| 扫描数据预处理 | `preprocess_laser_scan(raw_pts)` | Ch9 | ⭐⭐⭐ |
| 点云网格化 | `point_cloud_mesh(pts, resolution)` | Ch11预处理 | ⭐⭐⭐ |
| 三角剖分 | `triangulate_surface(p, q, U, V, PQw)` | Ch11 | ⭐⭐⭐ |

### Ch12 — 生产系统集成
| 算法 | 函数签名 | 依赖 | 难度 |
|------|---------|------|------|
| OpenGL渲染 | `render_curve_openGL(crv)` | Ch5 | ⭐⭐⭐ |
| DOOR格式导出 | `export_to_door_format(crv/surf)` | Ch12 | ⭐⭐ |
| MOTIF UI集成 | `motif_widget_create()` | Ch12 | ⭐⭐ |

## 头文件目录结构

```
include/nurbs/
├── core/
│   ├── concepts.hpp          # CurveConcept, SurfaceConcept
│   ├── types.hpp             # Point, Vector, KnotVector, WeightVector
│   ├── numeric.hpp           # PrecisionConfig, Tolerance
│   └── utilities.hpp         # index_helpers, span_queries
├── basis/
│   ├── bspline_basis.hpp     # Algorithm A3.1 (Ch4 core)
│   ├── knot_insertion.hpp    # Algorithm A3.2
│   ├── knot_refinement.hpp   # Algorithm A3.3
│   └── degree_elevation.hpp  # Algorithm A3.4
├── curve/
│   ├── nurbs_curve.hpp       # Ch5: NURBS curve class
│   ├── curve_derivatives.hpp # Algorithm A5.2
│   ├── curve_inversion.hpp   # Algorithm A5.6
│   └── curve_operations.hpp  # Ch7: subdivision, merging
├── surface/
│   ├── nurbs_surface.hpp     # Ch6: NURBS surface class
│   ├── surface_derivatives.hpp
│   ├── surface_insertion.hpp
│   └── surface_operations.hpp # Ch8: decomposition, smoothing
├── approximation/
│   ├── curve_approx.hpp      # Ch9
│   └── surface_approx.hpp    # Ch9
├── interpolation/
│   ├── global_interp.hpp     # Algorithm A10.1
│   ├── hermite_interp.hpp   # Algorithm A10.2
│   └── least_squares.hpp     # Algorithm A10.3/10.4
├── processing/
│   ├── laser_scan.hpp        # Ch11
│   └── point_cloud.hpp       # Ch11
├── visualization/
│   └── opengl_render.hpp     # Ch12
└── nurbs.hpp                 # Master include
```
