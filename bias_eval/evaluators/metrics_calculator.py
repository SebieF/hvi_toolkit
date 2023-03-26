import pandas as pd
import torch

from typing import List
from torchmetrics import Accuracy, Precision, Recall, F1Score, SpearmanCorrCoef, MatthewsCorrCoef, ConfusionMatrix

CLASSIFICATION_METRICS = sorted(["accuracy", "precision", "recall", "f1_score", "matthews-corr-coeff"])


def calculate_metric(predictions: List, targets: List, metric: str):
    num_classes = len(set(predictions))
    if metric == "accuracy":
        metric_function = Accuracy(task="binary", average="micro", num_classes=num_classes)
    elif metric == "precision":
        metric_function = Precision(task="binary", average="macro", num_classes=num_classes)
    elif metric == "recall":
        metric_function = Recall(task="binary", average="macro", num_classes=num_classes)
    elif metric == "f1_score":
        metric_function = F1Score(task="binary", average="macro", num_classes=num_classes)
    elif metric == "matthews-corr-coeff":
        metric_function = MatthewsCorrCoef(task="binary", num_classes=num_classes)
    elif metric == "confusion-matrix":
        metric_function = ConfusionMatrix(task="binary", num_classes=num_classes)
    else:
        raise NotImplementedError(f"Metric {metric} unknown!")
    return metric_function(torch.tensor(predictions).cpu().float(), torch.tensor(targets).cpu().float())


def calculate_all_metrics(predictions: List, targets: List, df_name: str) -> pd.DataFrame:
    metrics_df = pd.DataFrame(data={metric_key: [calculate_metric(predictions=predictions,
                                                                  targets=targets,
                                                                  metric=metric_key).item()]
                                    for metric_key in CLASSIFICATION_METRICS}).T
    metrics_df.name = df_name
    return metrics_df
