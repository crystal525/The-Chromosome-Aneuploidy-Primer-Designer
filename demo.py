#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
不平衡数据二分类模型演示
Demo for Imbalanced Binary Classification Model

这是一个简化的演示脚本，展示主要功能
"""

from imbalanced_binary_classifier import ImbalancedBinaryClassifier
import numpy as np

def main():
    print("=" * 60)
    print("不平衡数据二分类模型演示")
    print("=" * 60)
    
    # 创建分类器
    classifier = ImbalancedBinaryClassifier(random_state=42)
    
    # 生成不平衡数据
    print("\n1. 生成不平衡数据")
    print("-" * 30)
    X, y = classifier.generate_imbalanced_data(
        n_samples=2000,  # 较小的数据集用于快速演示
        n_features=15,
        class_ratio=0.08  # 8%的少数类
    )
    
    # 进行全面评估
    print("\n2. 进行全面评估")
    print("-" * 30)
    results_df, X_test, y_test = classifier.comprehensive_evaluation(X, y)
    
    # 显示最佳结果
    print("\n3. 最佳结果")
    print("-" * 30)
    best_result = results_df.loc[results_df['auc_score'].idxmax()]
    best_key = f"{best_result['model']}_{best_result['sampling_method']}"
    
    print(f"最佳模型组合: {best_key}")
    print(f"AUC Score: {best_result['auc_score']:.4f}")
    print(f"Average Precision: {best_result['ap_score']:.4f}")
    
    # 生成详细报告
    print("\n4. 详细分类报告")
    print("-" * 30)
    classifier.detailed_classification_report(X_test, y_test, best_key)
    
    print("\n" + "=" * 60)
    print("演示完成！")
    print("=" * 60)
    print("\n使用说明:")
    print("1. 运行完整程序: python3 imbalanced_binary_classifier.py")
    print("2. 运行使用示例: python3 example_usage.py")
    print("3. 查看详细文档: README_imbalanced_classifier.md")
    print("\n主要功能:")
    print("- 多种采样方法处理不平衡数据")
    print("- 多种分类算法对比")
    print("- 全面的评估指标")
    print("- 可视化结果展示")

if __name__ == "__main__":
    main()