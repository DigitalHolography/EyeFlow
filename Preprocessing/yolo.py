from ultralytics import YOLO

def import_model(model_path):
    """
    Imports and returns a YOLO model from the specified path.

    Args:
        model_path: path to a YOLO model (e.g., UltraLytics YOLOv11)

    Returns:
        YOLO model instance
    """

    model = YOLO(model_path)
    return model

def predict_segmentation(model, image, conf=0.25):
    """
    Predicts segmentation masks for an image.

    Args:
        model: the YOLO segmentation model
        image: input image for prediction
        conf: confidence threshold for predictions

    Returns:
        Predicted segmentation masks
    """

    results = model.predict(source=image, conf=conf, verbose=False, show=False, save=False)
    if not results or results[0].masks is None:
        print(f"[Warning] No mask found.")
        return None
    
    return results[0].masks

def predict_classification(model, image, conf=0.25):
    """
    Predicts classification labels and confidences for an image.

    Args:
        model: the YOLO classification model
        image: input image for prediction
        conf: confidence threshold for predictions

    Returns:
        Predicted classification labels and confidences
    """

    results = model.predict(source=image, conf=conf, verbose=False, show=False, save=False)
    if not results or results[0].probs is None:
        print(f"[Warning] No classification result found.")
        return None
    
    scores = results[0].probs.data.tolist()
    names = model.model.names

    # List of (label, confidence)
    predictions = [
        (names[i], scores[i])
        for i in range(len(scores))
        if scores[i] >= conf
    ]

    if not predictions:
        print(f"[Warning] No class above confidence threshold {conf}.")
        return None

    return predictions
