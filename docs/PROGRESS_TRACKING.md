# 进度看板 — NURBS Kernel

**更新时间**: 2026-06-03  
**总进度**: Ch4+Ch5 完成 (11/33 算法, 33%) | 剩余 22 算法

---

## 📋 看板视图

### ✅ 已完成 (Done)
| 章节 | 算法 | 文件 |
|------|------|------|
| Ch4 | bspline_basis | basis/bspline_basis.hpp |
| Ch4 | knot_insertion | basis/knot_insertion.hpp |
| Ch4 | knot_refinement | basis/knot_refinement.hpp |
| Ch4 | degree_elevation | basis/degree_elevation.hpp |
| Ch5 | nurbs_curve构造 | curve/nurbs_curve.hpp |
| Ch5 | curve_derivatives | curve/curve_derivatives.hpp |
| Ch5 | curve_inversion | curve/curve_inversion.hpp |

### 🎯 下一批次 — P0 核心 (可并行，4算法)
| 优先级 | 算法 | 依赖 | 难点 |
|--------|------|------|------|
| P0 | surface_derivatives | Ch6 | ⭐⭐⭐⭐ 偏导递归，2D数组返回 |
| P0 | surface_knot_insertion | Ch6 | ⭐⭐⭐ 方向参数处理 |
| P0 | construct_surface | Ch6 | ⭐⭐⭐ 张量积构造 |
| P0 | surface_degree_elevation | Ch6 | ⭐⭐⭐⭐ 双向升阶 |

### 📋 待开发 — P1 基础 (串行依赖P0)
| 章节 | 算法 | 依赖 | 难度 |
|------|------|------|------|
| Ch7 | curve_subdivide | Ch5 | ⭐⭐ |
| Ch7 | curve_merge | Ch7_subdivide | ⭐⭐ |
| Ch7 | reparametrize | Ch5 | ⭐⭐⭐ |
| Ch8 | surface_decompose | Ch6 | ⭐⭐⭐⭐ |
| Ch8 | surface_join | Ch8_decompose | ⭐⭐⭐⭐ |

### 📋 待开发 — P2 逼近 (依赖P1)
| 章节 | 算法 | 依赖 | 难度 |
|------|------|------|------|
| Ch7 | curve_fitting | Ch7 | ⭐⭐⭐ |
| Ch9 | curve_approximate | Ch7_fitting | ⭐⭐⭐ |
| Ch9 | surface_approximate | Ch9_curve | ⭐⭐⭐⭐ |
| Ch10 | global_interpolate | Ch5 | ⭐⭐⭐ |
| Ch10 | hermite_interpolate | Ch10_global | ⭐⭐⭐ |
| Ch10 | least_squares_approx | Ch9 | ⭐⭐⭐⭐ |

### 📋 待开发 — P3 应用 (依赖P2)
| 章节 | 算法 | 依赖 | 难度 |
|------|------|------|------|
| Ch11 | preprocess_laser_scan | Ch9 | ⭐⭐⭐ |
| Ch11 | point_cloud_mesh | Ch11_preprocess | ⭐⭐⭐ |
| Ch11 | triangulate_surface | Ch11 | ⭐⭐⭐ |
| Ch12 | render_curve_openGL | Ch5 | ⭐⭐⭐ |
| Ch12 | export_to_door_format | Ch12 | ⭐⭐ |
| Ch12 | motif_widget_create | Ch12 | ⭐⭐ |

---

## 🔗 依赖链分析

```
Ch5(完成)
  └─► Ch6construct_surface ─► Ch6derivatives ─► Ch6升阶
                              └─► Ch8decompose ─► Ch8光顺/拼接
                                              └─► Ch9曲线逼近 ─► Ch9曲面逼近
                                                           └─► Ch10全局插值 ─► Hermite插值 ─► 最小二乘
                                                                        └─► Ch11预处理 ─► 点云网格化 ─► 三角剖分
                                                                                        └─► Ch12渲染/导出/UI

Ch7拆分/合并/重参数化 (独立于Ch6链)
  └─► Ch7_curve_fitting ─► Ch9_curve_approx
```

**关键路径**: Ch6 → Ch8 → Ch9 → Ch10 → Ch11 → Ch12 (串行主干)

**可并行的分支**:
- Ch7曲线操作 (独立于Ch6链)
- Ch12 OpenGL导出 (仅依赖Ch5)

---

## ⚠️ 关键技术难点优先级

| 优先级 | 难点 | 原因 |
|--------|------|------|
| 🔴 P0 | **surface_derivatives 偏导递归** | Ch8/Ch10/Ch11全部依赖此算法，返回2D导数数组 |
| 🔴 P0 | **surface_degree_elevation 双向升阶** | 复杂度高，但解锁Ch8所有操作 |
| 🟡 P1 | **least_squares_approx 最小二乘** | 病态矩阵处理，迭代收敛 |
| 🟡 P1 | **surface_approximate 曲面逼近** | 需要曲线逼近 + 张量积 + 全局优化 |
| 🟡 P1 | **triangulate_surface 三角剖分** | NURBS转三角网格，精度控制 |
| 🟢 P2 | Ch11激光扫描噪声过滤 | 实际数据质量决定算法选择 |
| 🟢 P2 | OpenGL渲染器架构 | 与核心算法解耦，可后置 |

---

## 📁 代码质量要求

1. **Concept Checks**: 所有算法前加 `static_assert` 验证 CurveConcept/SurfaceConcept
2. **constexpr**: degree, span query 等可编译时计算的必须有 constexpr 版本
3. **编译时验证**: KnotVector 度数一致性、权重向量维度匹配
4. **测试覆盖**: 每算法至少对应一个 TEST_CASE，带边界条件
5. **精度容差**: 统一 Tolerance 类型，避免 magic number

---

## 📅 建议开发顺序

```
Phase 1 (当前): Ch6 曲面核心 (4算法, 2-3周)
  → construct_surface → surface_knot_insertion → surface_derivatives → surface_degree_elevation

Phase 2: Ch7曲线操作 + Ch8曲面分解 (7算法, 2-3周)
  → curve_subdivide → curve_merge → reparametrize → curve_fitting
  → surface_decompose → surface_smooth → surface_join

Phase 3: Ch9-Ch10 逼近与插值 (6算法, 3-4周)
  → curve_approximate → surface_approximate
  → global_interpolate → hermite_interpolate → least_squares_approx

Phase 4: Ch11-Ch12 应用集成 (6算法, 2-3周)
  → preprocess_laser_scan → point_cloud_mesh → triangulate_surface
  → render_openGL → export_door → motif_ui
```
