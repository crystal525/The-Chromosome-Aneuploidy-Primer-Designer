#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
不平衡数据二分类模型使用示例
Example Usage of Imbalanced Binary Classification Model

这个文件展示了如何使用不平衡二分类器的各种功能
"""

import numpy as np
import pandas as pd
from imbalanced_binary_classifier import ImbalancedBinaryClassifier

def example_1_basic_usage():
    """示例1: 基本使用方法"""
    print("=" * 60)
    print("示例1: 基本使用方法")
    print("=" * 60)
    
    # 创建分类器
    classifier = ImbalancedBinaryClassifier(random_state=42)
    
    # 生成不平衡数据
    X, y = classifier.generate_imbalanced_data(
        n_samples=5000,
        n_features=15,
        class_ratio=0.05  # 5%的少数类，非常不平衡
    )
    
    # 进行全面评估
    results_df, X_test, y_test = classifier.comprehensive_evaluation(X, y)
    
    # 找到最佳模型
    best_result = results_df.loc[results_df['auc_score'].idxmax()]
    best_key = f"{best_result['model']}_{best_result['sampling_method']}"
    
    print(f"\n最佳模型组合: {best_key}")
    print(f"AUC Score: {best_result['auc_score']:.4f}")
    print(f"Average Precision: {best_result['ap_score']:.4f}")


def example_2_custom_data():
    """示例2: 使用自定义数据"""
    print("\n" + "=" * 60)
    print("示例2: 使用自定义数据")
    print("=" * 60)
    
    # 创建分类器
    classifier = ImbalancedBinaryClassifier(random_state=42)
    
    # 模拟真实的不平衡数据
    np.random.seed(42)
    n_samples = 8000
    n_features = 25
    
    # 生成特征数据
    X = np.random.randn(n_samples, n_features)
    
    # 创建不平衡的标签（少数类占3%）
    minority_ratio = 0.03
    n_minority = int(n_samples * minority_ratio)
    n_majority = n_samples - n_minority
    
    y = np.concatenate([
        np.zeros(n_majority),  # 多数类
        np.ones(n_minority)    # 少数类
    ])
    
    # 打乱数据
    indices = np.random.permutation(n_samples)
    X = X[indices]
    y = y[indices]
    
    # 加载数据
    classifier.load_custom_data(X, y)
    
    # 评估
    results_df, X_test, y_test = classifier.comprehensive_evaluation(X, y)
    
    # 显示结果
    print("\n评估结果:")
    print(results_df.sort_values('auc_score', ascending=False).head())


def example_3_specific_sampling_methods():
    """示例3: 使用特定的采样方法"""
    print("\n" + "=" * 60)
    print("示例3: 使用特定的采样方法")
    print("=" * 60)
    
    classifier = ImbalancedBinaryClassifier(random_state=42)
    
    # 生成数据
    X, y = classifier.generate_imbalanced_data(
        n_samples=3000,
        class_ratio=0.08
    )
    
    # 划分数据
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    # 标准化
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # 只使用SMOTE和XGBoost
    from imblearn.over_sampling import SMOTE
    from sklearn.ensemble import RandomForestClassifier
    
    # SMOTE采样
    smote = SMOTE(random_state=42)
    X_smote, y_smote = smote.fit_resample(X_train_scaled, y_train)
    
    print(f"SMOTE后数据分布: {np.bincount(y_smote)}")
    
    # 训练XGBoost
    classifier.initialize_models()
    model, y_pred, y_pred_proba = classifier.train_and_evaluate(
        X_smote, X_test_scaled, y_smote, y_test,
        'SMOTE', 'XGBoost'
    )
    
    # 详细报告
    classifier.detailed_classification_report(
        X_test_scaled, y_test, 'XGBoost_SMOTE'
    )


def example_4_model_comparison():
    """示例4: 模型对比分析"""
    print("\n" + "=" * 60)
    print("示例4: 模型对比分析")
    print("=" * 60)
    
    classifier = ImbalancedBinaryClassifier(random_state=42)
    
    # 生成不同不平衡程度的数据
    imbalance_ratios = [0.1, 0.05, 0.02]  # 10%, 5%, 2%
    
    for ratio in imbalance_ratios:
        print(f"\n不平衡比例: {ratio:.1%}")
        print("-" * 30)
        
        X, y = classifier.generate_imbalanced_data(
            n_samples=4000,
            class_ratio=ratio
        )
        
        # 只评估几种关键组合
        from sklearn.model_selection import train_test_split
        from sklearn.preprocessing import StandardScaler
        from imblearn.over_sampling import SMOTE
        
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y
        )
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # 测试不同方法
        methods = [
            ('原始数据', X_train_scaled, y_train),
            ('SMOTE', *SMOTE(random_state=42).fit_resample(X_train_scaled, y_train))
        ]
        
        classifier.initialize_models()
        
        for method_name, X_method, y_method in methods:
            for model_name in ['逻辑回归', '随机森林', 'XGBoost']:
                try:
                    model, y_pred, y_pred_proba = classifier.train_and_evaluate(
                        X_method, X_test_scaled, y_method, y_test,
                        method_name, model_name
                    )
                except Exception as e:
                    print(f"  {model_name} + {method_name}: 失败")


def example_5_threshold_optimization():
    """示例5: 阈值优化"""
    print("\n" + "=" * 60)
    print("示例5: 阈值优化")
    print("=" * 60)
    
    classifier = ImbalancedBinaryClassifier(random_state=42)
    
    # 生成数据
    X, y = classifier.generate_imbalanced_data(
        n_samples=5000,
        class_ratio=0.06
    )
    
    # 训练最佳模型
    results_df, X_test, y_test = classifier.comprehensive_evaluation(X, y)
    best_result = results_df.loc[results_df['auc_score'].idxmax()]
    best_key = f"{best_result['model']}_{best_result['sampling_method']}"
    
    # 获取预测概率
    y_pred_proba = classifier.results[best_key]['y_pred_proba']
    
    # 尝试不同阈值
    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    
    print(f"\n最佳模型: {best_key}")
    print("阈值优化结果:")
    print("阈值    精确率    召回率    F1分数")
    print("-" * 35)
    
    for threshold in thresholds:
        y_pred_thresh = (y_pred_proba >= threshold).astype(int)
        
        from sklearn.metrics import precision_score, recall_score, f1_score
        
        precision = precision_score(y_test, y_pred_thresh, zero_division=0)
        recall = recall_score(y_test, y_pred_thresh, zero_division=0)
        f1 = f1_score(y_test, y_pred_thresh, zero_division=0)
        
        print(f"{threshold:.1f}     {precision:.4f}    {recall:.4f}    {f1:.4f}")


def main():
    """运行所有示例"""
    print("不平衡数据二分类模型使用示例")
    print("=" * 60)
    
    try:
        example_1_basic_usage()
        example_2_custom_data()
        example_3_specific_sampling_methods()
        example_4_model_comparison()
        example_5_threshold_optimization()
        
        print("\n" + "=" * 60)
        print("所有示例运行完成！")
        print("=" * 60)
        
    except Exception as e:
        print(f"运行示例时出错: {e}")
        print("请确保已安装所有依赖包:")
        print("pip install -r requirements.txt")


if __name__ == "__main__":
    main()