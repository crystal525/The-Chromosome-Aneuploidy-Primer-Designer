#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
不平衡数据二分类模型
Imbalanced Binary Classification Model

这个脚本提供了处理不平衡数据的完整解决方案，包括：
1. 数据生成和预处理
2. 多种处理不平衡数据的方法
3. 多种分类算法
4. 模型评估和可视化
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    classification_report, confusion_matrix, roc_auc_score, 
    precision_recall_curve, roc_curve, average_precision_score
)
from imblearn.over_sampling import SMOTE, ADASYN, BorderlineSMOTE
from imblearn.under_sampling import RandomUnderSampler, EditedNearestNeighbours
from imblearn.combine import SMOTEENN, SMOTETomek
import xgboost as xgb
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class ImbalancedBinaryClassifier:
    """不平衡数据二分类器"""
    
    def __init__(self, random_state=42):
        self.random_state = random_state
        self.scaler = StandardScaler()
        self.models = {}
        self.results = {}
        
    def generate_imbalanced_data(self, n_samples=10000, n_features=20, 
                                n_informative=None, n_redundant=None,
                                class_ratio=0.1, noise=0.1):
        """
        生成不平衡的二分类数据
        
        Parameters:
        -----------
        n_samples : int, 样本总数
        n_features : int, 特征数量
        n_informative : int, 有用特征数量
        n_redundant : int, 冗余特征数量
        class_ratio : float, 少数类比例
        noise : float, 噪声水平
        
        Returns:
        --------
        X, y : 特征矩阵和标签
        """
        print(f"生成不平衡数据...")
        print(f"样本总数: {n_samples}, 特征数: {n_features}")
        print(f"少数类比例: {class_ratio:.1%}")
        
        # 自动设置参数
        if n_informative is None:
            n_informative = min(15, n_features - 2)
        if n_redundant is None:
            n_redundant = min(5, n_features - n_informative - 1)
        
        # 确保参数合理
        n_informative = min(n_informative, n_features - 1)
        n_redundant = min(n_redundant, n_features - n_informative)
        
        print(f"有用特征: {n_informative}, 冗余特征: {n_redundant}")
        
        # 生成分类数据
        X, y = make_classification(
            n_samples=n_samples,
            n_features=n_features,
            n_informative=n_informative,
            n_redundant=n_redundant,
            n_clusters_per_class=1,
            weights=[1-class_ratio, class_ratio],  # 不平衡比例
            flip_y=noise,
            random_state=self.random_state
        )
        
        # 转换为DataFrame便于分析
        feature_names = [f'feature_{i}' for i in range(n_features)]
        self.X_df = pd.DataFrame(X, columns=feature_names)
        self.y_series = pd.Series(y, name='target')
        
        # 打印数据分布
        class_counts = np.bincount(y)
        print(f"类别分布: 多数类={class_counts[0]}, 少数类={class_counts[1]}")
        print(f"不平衡比例: {class_counts[1]/class_counts[0]:.3f}")
        
        return X, y
    
    def load_custom_data(self, X, y):
        """加载自定义数据"""
        self.X_df = pd.DataFrame(X)
        self.y_series = pd.Series(y, name='target')
        print(f"加载自定义数据: {X.shape[0]} 样本, {X.shape[1]} 特征")
        class_counts = np.bincount(y)
        print(f"类别分布: 多数类={class_counts[0]}, 少数类={class_counts[1]}")
        return X, y
    
    def apply_sampling_methods(self, X, y):
        """应用各种采样方法处理不平衡数据"""
        sampling_methods = {
            '原始数据': (X, y),
            'SMOTE': SMOTE(random_state=self.random_state),
            'ADASYN': ADASYN(random_state=self.random_state),
            'BorderlineSMOTE': BorderlineSMOTE(random_state=self.random_state),
            'RandomUnderSampler': RandomUnderSampler(random_state=self.random_state),
            'EditedNearestNeighbours': EditedNearestNeighbours(),
            'SMOTEENN': SMOTEENN(random_state=self.random_state),
            'SMOTETomek': SMOTETomek(random_state=self.random_state)
        }
        
        sampled_data = {}
        for name, method in sampling_methods.items():
            try:
                if name == '原始数据':
                    X_resampled, y_resampled = method
                else:
                    X_resampled, y_resampled = method.fit_resample(X, y)
                
                sampled_data[name] = (X_resampled, y_resampled)
                print(f"{name}: {X_resampled.shape[0]} 样本, "
                      f"类别分布 {np.bincount(y_resampled)}")
            except Exception as e:
                print(f"{name} 采样失败: {e}")
                
        return sampled_data
    
    def initialize_models(self):
        """初始化多种分类模型"""
        self.models = {
            '逻辑回归': LogisticRegression(
                random_state=self.random_state,
                max_iter=1000,
                class_weight='balanced'  # 自动平衡类别权重
            ),
            '随机森林': RandomForestClassifier(
                n_estimators=100,
                random_state=self.random_state,
                class_weight='balanced'
            ),
            'XGBoost': xgb.XGBClassifier(
                random_state=self.random_state,
                eval_metric='logloss',
                scale_pos_weight=1  # 将在训练时动态调整
            )
        }
        print("已初始化分类模型:")
        for name in self.models.keys():
            print(f"  - {name}")
    
    def train_and_evaluate(self, X_train, X_test, y_train, y_test, 
                          sampling_method='原始数据', model_name='逻辑回归'):
        """训练和评估模型"""
        print(f"\n训练模型: {model_name} + {sampling_method}")
        
        # 获取模型
        model = self.models[model_name].__class__(**self.models[model_name].get_params())
        
        # 对于XGBoost，计算类别权重
        if model_name == 'XGBoost':
            class_counts = np.bincount(y_train)
            scale_pos_weight = class_counts[0] / class_counts[1]
            model.set_params(scale_pos_weight=scale_pos_weight)
        
        # 训练模型
        model.fit(X_train, y_train)
        
        # 预测
        y_pred = model.predict(X_test)
        y_pred_proba = model.predict_proba(X_test)[:, 1]
        
        # 计算评估指标
        auc_score = roc_auc_score(y_test, y_pred_proba)
        ap_score = average_precision_score(y_test, y_pred_proba)
        
        # 存储结果
        key = f"{model_name}_{sampling_method}"
        self.results[key] = {
            'model': model,
            'y_pred': y_pred,
            'y_pred_proba': y_pred_proba,
            'auc_score': auc_score,
            'ap_score': ap_score,
            'confusion_matrix': confusion_matrix(y_test, y_pred)
        }
        
        print(f"AUC Score: {auc_score:.4f}")
        print(f"Average Precision: {ap_score:.4f}")
        
        return model, y_pred, y_pred_proba
    
    def comprehensive_evaluation(self, X, y, test_size=0.2):
        """全面评估所有采样方法和模型的组合"""
        print("=" * 60)
        print("开始全面评估...")
        
        # 划分训练测试集
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=self.random_state, 
            stratify=y
        )
        
        # 标准化特征
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # 应用采样方法
        sampled_data = self.apply_sampling_methods(X_train_scaled, y_train)
        
        # 初始化模型
        self.initialize_models()
        
        # 评估所有组合
        evaluation_results = []
        
        for sampling_name, (X_sampled, y_sampled) in sampled_data.items():
            for model_name in self.models.keys():
                try:
                    # 训练和评估
                    model, y_pred, y_pred_proba = self.train_and_evaluate(
                        X_sampled, X_test_scaled, y_sampled, y_test,
                        sampling_name, model_name
                    )
                    
                    # 记录结果
                    evaluation_results.append({
                        'sampling_method': sampling_name,
                        'model': model_name,
                        'auc_score': self.results[f"{model_name}_{sampling_name}"]['auc_score'],
                        'ap_score': self.results[f"{model_name}_{sampling_name}"]['ap_score']
                    })
                    
                except Exception as e:
                    print(f"评估失败 {model_name} + {sampling_name}: {e}")
        
        # 转换为DataFrame便于分析
        results_df = pd.DataFrame(evaluation_results)
        
        print("\n" + "=" * 60)
        print("评估结果汇总:")
        print("=" * 60)
        print(results_df.sort_values('auc_score', ascending=False).to_string(index=False))
        
        return results_df, X_test_scaled, y_test
    
    def plot_evaluation_results(self, results_df):
        """绘制评估结果"""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. AUC分数对比
        pivot_auc = results_df.pivot(index='sampling_method', columns='model', values='auc_score')
        sns.heatmap(pivot_auc, annot=True, fmt='.3f', cmap='viridis', ax=axes[0,0])
        axes[0,0].set_title('AUC Score 对比')
        axes[0,0].set_xlabel('模型')
        axes[0,0].set_ylabel('采样方法')
        
        # 2. Average Precision对比
        pivot_ap = results_df.pivot(index='sampling_method', columns='model', values='ap_score')
        sns.heatmap(pivot_ap, annot=True, fmt='.3f', cmap='plasma', ax=axes[0,1])
        axes[0,1].set_title('Average Precision 对比')
        axes[0,1].set_xlabel('模型')
        axes[0,1].set_ylabel('采样方法')
        
        # 3. 最佳模型组合的AUC分布
        best_combinations = results_df.nlargest(5, 'auc_score')
        axes[1,0].barh(range(len(best_combinations)), best_combinations['auc_score'])
        axes[1,0].set_yticks(range(len(best_combinations)))
        axes[1,0].set_yticklabels([f"{row['model']}\n{row['sampling_method']}" 
                                  for _, row in best_combinations.iterrows()])
        axes[1,0].set_xlabel('AUC Score')
        axes[1,0].set_title('Top 5 模型组合')
        
        # 4. 采样方法效果对比
        sampling_performance = results_df.groupby('sampling_method')['auc_score'].mean().sort_values(ascending=False)
        axes[1,1].bar(range(len(sampling_performance)), sampling_performance.values)
        axes[1,1].set_xticks(range(len(sampling_performance)))
        axes[1,1].set_xticklabels(sampling_performance.index, rotation=45, ha='right')
        axes[1,1].set_ylabel('平均 AUC Score')
        axes[1,1].set_title('采样方法平均效果')
        
        plt.tight_layout()
        plt.show()
    
    def plot_roc_curves(self, X_test, y_test, top_n=3):
        """绘制ROC曲线"""
        plt.figure(figsize=(10, 8))
        
        # 获取最佳的几个组合
        all_results = []
        for key, result in self.results.items():
            all_results.append({
                'key': key,
                'auc_score': result['auc_score']
            })
        
        best_results = sorted(all_results, key=lambda x: x['auc_score'], reverse=True)[:top_n]
        
        for result_info in best_results:
            key = result_info['key']
            result = self.results[key]
            
            fpr, tpr, _ = roc_curve(y_test, result['y_pred_proba'])
            auc_score = result['auc_score']
            
            plt.plot(fpr, tpr, label=f'{key} (AUC = {auc_score:.3f})')
        
        plt.plot([0, 1], [0, 1], 'k--', label='随机分类器')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('假正率 (False Positive Rate)')
        plt.ylabel('真正率 (True Positive Rate)')
        plt.title('ROC 曲线对比')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()
    
    def plot_precision_recall_curves(self, X_test, y_test, top_n=3):
        """绘制精确率-召回率曲线"""
        plt.figure(figsize=(10, 8))
        
        # 获取最佳的几个组合
        all_results = []
        for key, result in self.results.items():
            all_results.append({
                'key': key,
                'ap_score': result['ap_score']
            })
        
        best_results = sorted(all_results, key=lambda x: x['ap_score'], reverse=True)[:top_n]
        
        for result_info in best_results:
            key = result_info['key']
            result = self.results[key]
            
            precision, recall, _ = precision_recall_curve(y_test, result['y_pred_proba'])
            ap_score = result['ap_score']
            
            plt.plot(recall, precision, label=f'{key} (AP = {ap_score:.3f})')
        
        # 计算随机分类器的AP
        random_ap = np.sum(y_test) / len(y_test)
        plt.axhline(y=random_ap, color='k', linestyle='--', label=f'随机分类器 (AP = {random_ap:.3f})')
        
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('召回率 (Recall)')
        plt.ylabel('精确率 (Precision)')
        plt.title('精确率-召回率曲线对比')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()
    
    def detailed_classification_report(self, X_test, y_test, best_key):
        """生成详细的分类报告"""
        if best_key not in self.results:
            print(f"未找到结果: {best_key}")
            return
        
        result = self.results[best_key]
        y_pred = result['y_pred']
        
        print(f"\n最佳模型组合: {best_key}")
        print("=" * 50)
        print("分类报告:")
        print(classification_report(y_test, y_pred, target_names=['多数类', '少数类']))
        
        print("\n混淆矩阵:")
        cm = result['confusion_matrix']
        print(f"                预测")
        print(f"实际    多数类  少数类")
        print(f"多数类    {cm[0,0]:4d}   {cm[0,1]:4d}")
        print(f"少数类    {cm[1,0]:4d}   {cm[1,1]:4d}")
        
        # 计算关键指标
        tn, fp, fn, tp = cm.ravel()
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        
        print(f"\n关键指标:")
        print(f"精确率 (Precision): {precision:.4f}")
        print(f"召回率 (Recall): {recall:.4f}")
        print(f"F1分数: {f1:.4f}")
        print(f"特异性 (Specificity): {specificity:.4f}")
        print(f"AUC Score: {result['auc_score']:.4f}")
        print(f"Average Precision: {result['ap_score']:.4f}")


def main():
    """主函数 - 演示如何使用不平衡二分类器"""
    print("不平衡数据二分类模型演示")
    print("=" * 50)
    
    # 创建分类器实例
    classifier = ImbalancedBinaryClassifier(random_state=42)
    
    # 生成不平衡数据
    X, y = classifier.generate_imbalanced_data(
        n_samples=10000,
        n_features=20,
        class_ratio=0.1  # 10%的少数类
    )
    
    # 进行全面评估
    results_df, X_test, y_test = classifier.comprehensive_evaluation(X, y)
    
    # 绘制评估结果
    classifier.plot_evaluation_results(results_df)
    
    # 绘制ROC曲线
    classifier.plot_roc_curves(X_test, y_test, top_n=3)
    
    # 绘制精确率-召回率曲线
    classifier.plot_precision_recall_curves(X_test, y_test, top_n=3)
    
    # 找到最佳模型组合
    best_result = results_df.loc[results_df['auc_score'].idxmax()]
    best_key = f"{best_result['model']}_{best_result['sampling_method']}"
    
    # 生成详细报告
    classifier.detailed_classification_report(X_test, y_test, best_key)
    
    print("\n" + "=" * 50)
    print("演示完成！")
    print("你可以使用以下方法自定义数据:")
    print("1. classifier.load_custom_data(X, y) - 加载自己的数据")
    print("2. classifier.comprehensive_evaluation(X, y) - 进行全面评估")
    print("3. 查看 results_df 获取所有结果")


if __name__ == "__main__":
    main()