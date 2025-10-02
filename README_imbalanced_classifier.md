# 不平衡数据二分类模型

这是一个完整的Python解决方案，用于处理不平衡数据的二分类问题。该模型集成了多种采样方法、分类算法和评估指标，帮助你在不平衡数据集上构建高性能的分类器。

## 功能特点

### 🔧 数据处理
- **数据生成**: 自动生成不平衡的二分类数据
- **自定义数据**: 支持加载你自己的数据集
- **数据标准化**: 自动进行特征标准化处理

### 📊 采样方法
- **过采样**: SMOTE, ADASYN, BorderlineSMOTE
- **欠采样**: RandomUnderSampler, EditedNearestNeighbours
- **混合采样**: SMOTEENN, SMOTETomek
- **权重调整**: 自动类别权重平衡

### 🤖 分类算法
- **逻辑回归**: 带类别权重平衡
- **随机森林**: 集成学习方法
- **XGBoost**: 梯度提升算法

### 📈 评估指标
- **AUC Score**: ROC曲线下面积
- **Average Precision**: 精确率-召回率曲线下面积
- **混淆矩阵**: 详细的分类结果
- **精确率、召回率、F1分数**: 关键性能指标

### 📊 可视化
- **热力图**: 不同方法组合的效果对比
- **ROC曲线**: 多个模型的ROC曲线对比
- **精确率-召回率曲线**: PR曲线分析
- **性能排名**: 最佳模型组合展示

## 安装依赖

```bash
pip install -r requirements.txt
```

## 快速开始

### 基本使用

```python
from imbalanced_binary_classifier import ImbalancedBinaryClassifier

# 创建分类器
classifier = ImbalancedBinaryClassifier(random_state=42)

# 生成不平衡数据
X, y = classifier.generate_imbalanced_data(
    n_samples=10000,
    n_features=20,
    class_ratio=0.1  # 10%的少数类
)

# 进行全面评估
results_df, X_test, y_test = classifier.comprehensive_evaluation(X, y)

# 查看结果
print(results_df.sort_values('auc_score', ascending=False))
```

### 使用自定义数据

```python
import numpy as np
import pandas as pd

# 准备你的数据
X = your_features  # 特征矩阵
y = your_labels    # 标签 (0 和 1)

# 加载数据
classifier = ImbalancedBinaryClassifier()
classifier.load_custom_data(X, y)

# 评估
results_df, X_test, y_test = classifier.comprehensive_evaluation(X, y)
```

### 可视化结果

```python
# 绘制评估结果热力图
classifier.plot_evaluation_results(results_df)

# 绘制ROC曲线
classifier.plot_roc_curves(X_test, y_test, top_n=3)

# 绘制精确率-召回率曲线
classifier.plot_precision_recall_curves(X_test, y_test, top_n=3)
```

## 详细使用指南

### 1. 数据生成参数

```python
X, y = classifier.generate_imbalanced_data(
    n_samples=10000,      # 样本总数
    n_features=20,        # 特征数量
    n_informative=15,     # 有用特征数量
    n_redundant=5,        # 冗余特征数量
    class_ratio=0.1,      # 少数类比例 (0.1 = 10%)
    noise=0.1             # 噪声水平
)
```

### 2. 采样方法说明

| 方法 | 类型 | 描述 |
|------|------|------|
| SMOTE | 过采样 | 合成少数类过采样技术 |
| ADASYN | 过采样 | 自适应合成采样 |
| BorderlineSMOTE | 过采样 | 边界SMOTE |
| RandomUnderSampler | 欠采样 | 随机欠采样 |
| EditedNearestNeighbours | 欠采样 | 编辑最近邻 |
| SMOTEENN | 混合 | SMOTE + 编辑最近邻 |
| SMOTETomek | 混合 | SMOTE + Tomek链接 |

### 3. 模型参数调整

```python
# 自定义模型参数
classifier.models = {
    '逻辑回归': LogisticRegression(
        random_state=42,
        max_iter=1000,
        class_weight='balanced',
        C=1.0  # 正则化参数
    ),
    '随机森林': RandomForestClassifier(
        n_estimators=200,  # 树的数量
        max_depth=10,      # 最大深度
        random_state=42,
        class_weight='balanced'
    ),
    'XGBoost': xgb.XGBClassifier(
        n_estimators=100,
        max_depth=6,
        learning_rate=0.1,
        random_state=42
    )
}
```

### 4. 阈值优化

```python
# 获取最佳模型的预测概率
best_key = "XGBoost_SMOTE"
y_pred_proba = classifier.results[best_key]['y_pred_proba']

# 尝试不同阈值
thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
for threshold in thresholds:
    y_pred_thresh = (y_pred_proba >= threshold).astype(int)
    # 计算性能指标...
```

## 评估指标解释

### AUC Score
- **范围**: 0-1
- **含义**: ROC曲线下面积，衡量模型区分正负样本的能力
- **解释**: 
  - 0.5: 随机分类器
  - 0.7-0.8: 较好
  - 0.8-0.9: 很好
  - >0.9: 优秀

### Average Precision
- **范围**: 0-1
- **含义**: 精确率-召回率曲线下面积
- **特点**: 对不平衡数据更敏感，更适合评估不平衡分类

### 混淆矩阵
```
                预测
实际    多数类  少数类
多数类    TN     FP
少数类    FN     TP
```

- **TN**: 真负例 (正确预测为多数类)
- **FP**: 假正例 (错误预测为少数类)
- **FN**: 假负例 (错误预测为多数类)
- **TP**: 真正例 (正确预测为少数类)

## 最佳实践建议

### 1. 数据不平衡程度
- **轻度不平衡** (1:3 到 1:10): 使用类别权重调整
- **中度不平衡** (1:10 到 1:100): 使用SMOTE等过采样
- **重度不平衡** (>1:100): 使用混合采样方法

### 2. 模型选择
- **线性关系**: 逻辑回归 + 权重调整
- **非线性关系**: 随机森林或XGBoost
- **大数据集**: XGBoost通常表现更好

### 3. 评估策略
- **主要指标**: AUC Score 和 Average Precision
- **业务指标**: 根据业务需求关注精确率或召回率
- **阈值优化**: 根据业务成本调整分类阈值

### 4. 交叉验证
```python
from sklearn.model_selection import StratifiedKFold

# 使用分层K折交叉验证
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = cross_val_score(model, X, y, cv=skf, scoring='roc_auc')
print(f"交叉验证AUC: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
```

## 常见问题

### Q: 如何处理极度不平衡的数据？
A: 对于极度不平衡的数据（如1:1000），建议：
1. 使用SMOTEENN或SMOTETomek混合采样
2. 考虑异常检测方法
3. 使用集成学习方法

### Q: 如何选择合适的采样方法？
A: 建议按以下顺序尝试：
1. 原始数据 + 类别权重
2. SMOTE
3. BorderlineSMOTE
4. SMOTEENN
5. 根据结果选择最佳方法

### Q: 如何解释模型结果？
A: 重点关注：
1. AUC Score > 0.7 通常是可以接受的
2. Average Precision 对不平衡数据更敏感
3. 混淆矩阵帮助理解错误类型
4. ROC和PR曲线显示不同阈值下的性能

## 运行示例

```bash
# 运行主程序
python imbalanced_binary_classifier.py

# 运行使用示例
python example_usage.py
```

## 文件结构

```
workspace/
├── imbalanced_binary_classifier.py  # 主程序
├── example_usage.py                 # 使用示例
├── requirements.txt                 # 依赖包
└── README_imbalanced_classifier.md  # 说明文档
```

## 许可证

MIT License

## 贡献

欢迎提交Issue和Pull Request来改进这个项目！