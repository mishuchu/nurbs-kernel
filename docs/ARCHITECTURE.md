# NURBS Kernel 架构设计文档

> 版本：v0.1 | 日期：2026-06-03 | 状态：草稿

## 1. 项目概述

- **目标**：复现 The NURBS Book (Piegl & Tiller, 2nd ed.) 第4-11章算法
- **风格**：C++20 header-only，纯模板，无依赖
- **当前状态**：Ch4(4/5), Ch5(3/6), Ch6-Ch12 未实现

---

## 2. 代码审查（已完成部分）

### 2.1 bspline_basis.hpp

**问题清单：**

| 严重度 | 位置 | 问题 |
|--------|------|------|
| 🔴 高 | L151-155 | `compute_basis_function_derivatives` 是 stub — 导数表计算为空循环，未实现 Algorithm A3.2 |
| 🟡 中 | L39-42 | `Tolerance<T>` 使用 `defaults()` 但传参方式不一致 — 有调用传，有不传 |
| 🟡 中 | L68-83 | `basis_function` 递归实现 Cox-de Boor，无 memoization，深度 p 可能导致栈溢出 |
| 🟡 中 | L92-96 | `BasisFunctionDerivatives` 注释写 Algorithm A3.2，但实际代码未实现导数计算 |
| 🟢 低 | L35 | 返回 `BasisFunctionValues<T>` 但实际 span 索引由内部 compute 返回，未显式使用 |

**结论**：基础求值正确，导数计算未完成。

---

### 2.2 knot_insertion.hpp

**问题清单：**

| 严重度 | 位置 | 问题 |
|--------|------|------|
| 🟡 中 | L68-69 | 多重性检查 `U[i] == u` 用精确比较，缺少容差 — 浮点场景可能漏判 |
| 🟡 中 | L81-86 | 插入位置逻辑：注释说"after position k"但循环是 `<= k` 后插，算法正确但注释误导 |
| 🟡 中 | L98-102 | alpha 计算中若分母为 0 未做安全处理（理论上不会发生但防御性编程缺失）|
| 🟢 低 | L146-153 | `curve_knot_insertion` 是 `insert_knot` 的 thin wrapper，建议合并或明确区分用途 |

**结论**：算法核心正确，容差处理缺失。

---

### 2.3 degree_elevation.hpp

**问题清单：**

| 严重度 | 位置 | 问题 |
|--------|------|------|
| 🔴 高 | L95-118 | `elevate_degree` 未实现 — 只构建了 knot vector 结构，`Qw_bar` 全零，控制点未计算 |
| 🔴 高 | L77-82 | 二项式系数 binom[i][j]=1 是占位符，实际计算依赖 bezier_elevate_degree 中的递归公式 |
| 🟡 中 | L88-93 | 注释说"simplified"和"placeholder"，但无 TODO/FIXME 标记，生产代码易忽略 |
| 🟡 中 | L84-91 | 循环逻辑在 `p_bar` 变化时控制点数量不匹配 |

**结论**：`bezier_elevate_degree` 正确，但主函数 `elevate_degree` 未完成。

---

### 2.4 knot_refinement.hpp

**问题清单：**

| 严重度 | 位置 | 问题 |
|--------|------|------|
| 🔴 高 | L86-92 | `refine_knot_vector` 控制点数组全零占位，未从 `Qw` 计算 |
| 🟡 中 | L101-156 | `refine_knot_vector_full` 有完整算法但实现有 bug — L144 循环范围计算错误 |
| 🟡 中 | L119 | `find_span` 用 `X.front()` 但未验证 X 非空（r==0 时在 L60 处理了，但多了一步）|

**结论**：`refine_knot_vector_full` 接近正确但有索引 bug。

---

### 2.5 curve_derivatives.hpp

**问题清单：**

| 严重度 | 位置 | 问题 |
|--------|------|------|
| 🟡 中 | L93-103 | 导数公式索引可能越界：`U[k + p + ri - i + 1]` 当 i>p 时访问越界 |
| 🟢 低 | L125-130 | 重载 `curve_derivatives(u, U, p, Pw)` 调用自身 d=p，默认计算所有导数，逻辑正确 |

**结论**：主体正确，边界索引需审查。

---

### 2.6 curve_inversion.hpp

**问题清单：**

| 严重度 | 位置 | 问题 |
|--------|------|------|
| 🔴 高 | L87-88 | 硬编码 2D 运算：`dot_diff_cp` 和 `denom` 只取 [0] 和 [1] 坐标，不支持 3D |
| 🟡 中 | L64-66 | `is_valid_knot_vector` 返回 pair，但后续用 `U[p]` 和 `U[U.size()-1-p]` 不安全 |
| 🟡 中 | L146-176 | `find_initial_guess_candidates` 只采样 20 个点且不考虑 knot spans，精度不足 |
| 🟡 中 | L117-122 | `curve_inversion_auto_guess` 用中点作为初始猜测，对弯折曲线可能不收敛 |

**结论**：仅 2D 支持，3D 扩展需重写。

---

### 2.7 nurbs_curve.hpp

**问题清单：**

| 严重度 | 位置 | 问题 |
|--------|------|------|
| 🔴 高 | L53-58 | `CurveConcept` 要求 `weights()` 方法，但 `NURBSCurve` 未暴露公开接口 |
| 🟡 中 | L162-193 | `evaluate_derivatives` 实现是简化版，非 Algorithm A5.2 — 无控制点链式计算 |
| 🟡 中 | L109 | `if (u == u_max)` 精确比较，浮点边界风险 |
| 🟢 低 | L225-240 | `validate()` 检查 `is_valid_knot_vector` 但此函数要求 clamped 曲线，非通用 |

**结论**：核心类可用但 derivative 实现需按算法重构，concept 一致性需修复。

---

### 2.8 concepts.hpp

**问题清单：**

| 严重度 | 位置 | 问题 |
|--------|------|------|
| 🔴 高 | L53-58 | `CurveConcept` 需要 `weights()` 但 `NURBSCurve` 未实现 — concept 违背 |
| 🟡 中 | L48-58 | `CurveConcept` 中 `weights()` 返回 `WeightVector<T>` 但无 `NURBSPoint::weight()` 接口 |
| 🟡 中 | L66-81 | `SurfaceConcept` 定义完整但 `NURBSSurface` 只有 stub — 无法编译期检查 |

**结论**：concept 定义与实现不匹配，需同步修复。

---

## 3. API 一致性分析

### 3.1 函数签名不一致

| 模块 | 签名模式 |
|------|----------|
| `basis/` | `fn(n, p, u, U, ...)` — n=控制点数-1，显式传参 |
| `curve/` | `fn(u, p, U, Pw, ...)` — 从 Pw 隐式算 n |
| `surface/` | 未实现，待定义 |

**建议**：统一为 `fn(U, p, Pw, ...)` — 去掉显式 n，由 `Pw.size()-1` 推导。

### 3.2 命名不一致

| 当前 | 建议 | 说明 |
|------|------|------|
| `BasisFunctionValues` | `BasisFunctions` | 简洁 |
| `KnotInsertionResult` | `KnotInsertionResult` | OK |
| `CurveKnotInsertionResult` | 合并到 `KnotInsertionResult` | 重复 |
| `DegreeElevationResult` | `DegreeElevationResult` | OK |

### 3.3 Result 类型命名约定

所有 Result struct 统一为：
```cpp
template <NumericScalar_ T>
struct <Operation>Result {
    // 输出字段
};
```

---

## 4. C++20 特性使用情况

| 特性 | 使用 | 评价 |
|------|------|------|
| `concept` | `NumericScalar_`, `ParametricEntity`, `CurveConcept`, `SurfaceConcept` | ✅ 正确使用 |
| `requires` | 仅在 concept 定义中 | ✅ 简洁 |
| `[[nodiscard]]` | 大部分函数 | ✅ 一致 |
| `constexpr` | `Tolerance`, `Point`, `PrecisionConfig` 部分 | ⚠️ 不完整，核心算法未 constexpr |
| `static constexpr` | `Tolerance::defaults()`, `PrecisionConfig` | ✅ |

---

## 5. Ch6-Ch12 目录结构与 API 设计

### 5.1 目录结构

```
include/nurbs/
├── core/
│   ├── concepts.hpp          # CurveConcept, SurfaceConcept (需修复)
│   ├── types.hpp             # Point, NURBSPoint, KnotVector, WeightVector
│   ├── numeric.hpp           # PrecisionConfig, Tolerance, NumericScalar_
│   └── utilities.hpp         # find_span, compute_basis_functions
├── basis/
│   ├── bspline_basis.hpp     # A3.1 ✅(导数stub)
│   ├── knot_insertion.hpp    # A3.2 ✅
│   ├── knot_refinement.hpp   # A3.3 ⚠️(有bug)
│   ├── degree_elevation.hpp  # A3.4 🔴(未完成)
│   └── knot_removal.hpp      # A3.5 ❌(未实现)
├── curve/
│   ├── nurbs_curve.hpp       # A5.1 ✅
│   ├── curve_derivatives.hpp # A5.2 ✅(需审查边界)
│   ├── curve_inversion.hpp   # A5.6 🔴(仅2D)
│   ├── curve_knot_insertion.hpp # A5.1 wrapper ❌(未实现)
│   ├── curve_degree_elevation.hpp # A5.3 ❌(未实现)
│   └── curve_integrate.hpp  # ❌(未实现)
├── surface/
│   ├── nurbs_surface.hpp     # stub → A6.1
│   ├── surface_derivatives.hpp # A6.2 ❌(未实现)
│   ├── surface_knot_insertion.hpp # A6.3 ❌(未实现)
│   └── surface_degree_elevation.hpp # A6.4 ❌(未实现)
├── operations/                # Ch7: 曲线操作
│   ├── curve_subdivide.hpp   # A7.1 ❌
│   ├── curve_merge.hpp       # A7.2 ❌
│   ├── curve_reparametrize.hpp # A7.3 ❌
│   └── curve_fitting.hpp     # A7.4 ❌
├── approximation/            # Ch9: 逼近
│   ├── curve_approx.hpp      # A9.1/A9.2 ❌
│   └── surface_approx.hpp    # A9.3 ❌
├── interpolation/            # Ch10: 插值
│   ├── global_interp.hpp     # A10.1 ❌
│   ├── hermite_interp.hpp    # A10.2 ❌
│   └── least_squares.hpp     # A10.3/10.4 ❌
├── processing/              # Ch11: 数据处理
│   ├── laser_scan.hpp        # ❌
│   └── point_cloud.hpp       # ❌
├── visualization/           # Ch12: 生产系统
│   └── opengl_render.hpp    # ❌
└── nurbs.hpp               # Master include
```

### 5.2 API 设计原则

#### 原则 1：签名统一

```cpp
// 基础函数
template <NumericScalar_ T>
[[nodiscard]] ResultType<T> operation(const KnotVector<T>& U, int p,
                                       const std::vector<NURBSPoint<T>>& Pw, Args...);

// 替代当前的不一致签名
// ❌ 当前: insert_knot(T u, int p, const KnotVector<T>& U, const std::vector<NURBSPoint<T>>& Pw)
// ✅ 建议: insert_knot(const KnotVector<T>& U, int p, T u, const std::vector<NURBSPoint<T>>& Pw)
```

#### 原则 2：Result 类型明确

```cpp
template <NumericScalar_ T>
struct KnotInsertionResult {
    KnotVector<T> knot_vector;           // U_bar
    std::vector<NURBSPoint<T>> control_points; // Qw_bar
    std::size_t span;                     // k
};

// 曲面结果增加方向参数
template <NumericScalar_ T>
struct SurfaceKnotInsertionResult {
    KnotVector<T> u_knot_vector;
    KnotVector<T> v_knot_vector;
    std::vector<std::vector<NURBSPoint<T>>> control_points;
    char direction; // 'u' or 'v'
};
```

#### 原则 3：使用 requires 约束

```cpp
template <typename T>
concept CurveEvaluator = requires(T c, double u) {
    { c.evaluate(u) } -> std::same_as<NURBSPoint<typename T::scalar_type>>;
    { c.degree() } -> std::same_as<int>;
    { c.knot_vector() } -> std::same_as<const KnotVector<typename T::scalar_type>&>;
};

template <CurveEvaluator C, NumericScalar_ T>
[[nodiscard]] T curve_inversion(const Point<T>& target, const C& curve, Tolerance<T> tol = {});
```

#### 原则 4：容差管理

```cpp
// 显式容差参数，无默认值时使用 Tolerance<T>::defaults()
// 避免隐式默认值导致的不同步
template <NumericScalar_ T>
[[nodiscard]] bool knots_equal(T a, T b, Tolerance<T> tol) {
    return tol.eq(a, b);
}
```

#### 原则 5：算法标注

```cpp
// 文件头注释格式
// algorithm_name.hpp — Algorithm A?.? (NURBS Book, 2nd ed., p.???)
// 中文功能描述
#pragma once

namespace nurbs::<module> {
    // 算法实现
}
```

---

## 6. Ch6-Ch12 API 规范

### 6.1 Ch6 — NURBS 曲面

```cpp
namespace nurbs::surface {

// A6.1 — 曲面构造
template <NumericScalar_ T>
class NURBSSurface {
public:
    using scalar_type = T;

    NURBSSurface() = default;
    NURBSSurface(int p, int q, KnotVector<T> U, KnotVector<T> V,
                 std::vector<std::vector<NURBSPoint<T>>> PQw);

    [[nodiscard]] int u_degree() const noexcept;
    [[nodiscard]] int v_degree() const noexcept;
    [[nodiscard]] const KnotVector<T>& u_knot_vector() const noexcept;
    [[nodiscard]] const KnotVector<T>& v_knot_vector() const noexcept;
    [[nodiscard]] const std::vector<std::vector<NURBSPoint<T>>>& control_points() const noexcept;
    [[nodiscard]] std::pair<T,T> parameter_domain_u() const noexcept;
    [[nodiscard]] std::pair<T,T> parameter_domain_v() const noexcept;

    [[nodiscard]] NURBSPoint<T> evaluate(T u, T v) const;
    [[nodiscard]] std::vector<NURBSPoint<T>> evaluate_derivatives(T u, T v, int d) const;

    // A6.3 — 节点插入
    [[nodiscard]] NURBSSurface insert_knot_u(T u) const;
    [[nodiscard]] NURBSSurface insert_knot_v(T v) const;

    // A6.4 — 升阶
    [[nodiscard]] NURBSSurface elevate_u_degree(int t) const;
    [[nodiscard]] NURBSSurface elevate_v_degree(int t) const;
};

// A6.2 — 曲面导数
template <NumericScalar_ T>
struct SurfaceDerivatives {
    int u_degree, v_degree;
    std::vector<std::vector<std::vector<NURBSPoint<T>>>> derivatives; // [du][dv][point]
};

template <NumericScalar_ T>
[[nodiscard]] SurfaceDerivatives<T>
surface_derivatives(T u, T v, int p, int q,
                    const KnotVector<T>& U, const KnotVector<T>& V,
                    const std::vector<std::vector<NURBSPoint<T>>>& PQw,
                    int du, int dv);

} // namespace nurbs::surface
```

### 6.2 Ch7 — 曲线操作

```cpp
namespace nurbs::operations {

// A7.1 — 曲线拆分
template <NumericScalar_ T>
struct CurveSubdivisionResult {
    NURBSCurve<T> left;   // [u0, u]
    NURBSCurve<T> right;  // [u, u1]
    T split_point;
};

template <NumericScalar_ T>
[[nodiscard]] CurveSubdivisionResult<T>
curve_subdivide(T u, const NURBSCurve<T>& curve);

// A7.2 — 曲线合并
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
curve_merge(const NURBSCurve<T>& c1, const NURBSCurve<T>& c2);

// A7.3 — 重参数化
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
curve_reparametrize(const KnotVector<T>& U_new, const NURBSCurve<T>& curve);

// A7.4 — 曲线拟合
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
curve_fitting(const std::vector<Point<T>>& points, int p, int num_control_points);

} // namespace nurbs::operations
```

### 6.3 Ch9 — 逼近

```cpp
namespace nurbs::approximation {

// A9.1/A9.2 — 曲线逼次
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
curve_approximate(const std::vector<Point<T>>& points, int p, int num_knots);

// A9.3 — 曲面逼次
template <NumericScalar_ T>
[[nodiscard]] NURBSSurface<T>
surface_approximate(const std::vector<std::vector<Point<T>>>& points,
                    int u_degree, int v_degree,
                    const KnotVector<T>& U, const KnotVector<T>& V);

} // namespace nurbs::approximation
```

### 6.4 Ch10 — 插值

```cpp
namespace nurbs::interpolation {

// A10.1 — 全局插值
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
global_interpolate(const std::vector<Point<T>>& points, int p);

// A10.2 — Hermite 插值
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
hermite_interpolate(const Point<T>& P0, const Point<T>& T0,
                    const Point<T>& P1, const Point<T>& T1,
                    int p);

// A10.3/A10.4 — 最小二乘
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
least_squares_approx(const std::vector<Point<T>>& points, int p, int n);

} // namespace nurbs::interpolation
```

### 6.5 Ch11 — 数据处理

```cpp
namespace nurbs::processing {

// 激光扫描预处理
template <NumericScalar_ T>
[[nodiscard]] std::vector<Point<T>>
preprocess_laser_scan(const std::vector<Point<T>>& raw, Tolerance<T> tol);

// 点云网格化
template <NumericScalar_ T>
[[nodiscard]] std::vector<std::vector<Point<T>>>
point_cloud_mesh(const std::vector<Point<T>>& cloud, T resolution);

// 曲面的三角剖分
template <NumericScalar_ T>
[[nodiscard]] std::vector<std::array<std::size_t, 3>>
triangulate(const NURBSSurface<T>& surf, std::size_t u_segments, std::size_t v_segments);

} // namespace nurbs::processing
```

### 6.6 Ch12 — 生产系统

```cpp
namespace nurbs::visualization {

// OpenGL 渲染辅助
template <NumericScalar_ T>
void render_curve(const NURBSCurve<T>& curve, std::size_t num_samples);

// DOOR 格式导出
template <NumericScalar_ T>
[[nodiscard]] std::string export_to_door(const NURBSCurve<T>& curve);

template <NumericScalar_ T>
[[nodiscard]] std::string export_to_door(const NURBSSurface<T>& surf);

} // namespace nurbs::visualization
```

---

## 7. 修复优先级

| 优先级 | 文件 | 问题 | 修复工作量 |
|--------|------|------|------------|
| P0 | concepts.hpp | `CurveConcept` 缺少 `weights()` | 1h |
| P0 | curve_inversion.hpp | 仅 2D 支持 | 2h |
| P1 | degree_elevation.hpp | `elevate_degree` 未实现 | 4h |
| P1 | bspline_basis.hpp | `compute_basis_function_derivatives` stub | 3h |
| P1 | knot_refinement.hpp | `refine_knot_vector_full` 索引 bug | 2h |
| P2 | nurbs_curve.hpp | `evaluate_derivatives` 简化实现 | 3h |
| P2 | curve_derivatives.hpp | 导数边界索引 | 1h |

---

## 8. 总结

**已完成部分质量**：基础结构良好，核心算法大部分正确，但存在多处 stub/placeholder 和边界 case 处理缺失。

**主要技术债**：
1. `degree_elevation.hpp` 主函数未实现
2. `curve_inversion.hpp` 仅 2D，扩展需重写
3. Concept 定义与实现不匹配
4. 函数签名不一致（n vs Pw.size()-1）

**下一步**：
1. 修复 P0 问题（concept 和 inversion）使代码可编译
2. 完成 degree_elevation 和 basis derivatives
3. 按第6节 API 规范实现 Ch6-Ch12

---

*本文档由 Architect Agent 自动生成，审查了 7 个头文件，发现 P0 问题 2 个，P1 问题 3 个，P2 问题 2 个。*