def yolo_predict(model, image, conf=0.25):
    """
    Predicts segmentation masks for multiple images and displays
    the exact predicted mask contours (borders) for each image.

    Args:
        model_path: path to the YOLO segmentation model (e.g., YOLOv8)
        image: input image for prediction
        conf: confidence threshold for predictions
    """

    # Run prediction
    results = model.predict(source=image, conf=conf, verbose=False, show=False, save=False)
    if not results or results[0].masks is None:
        print(f"[Warning] No mask found.")
        return None
    return results[0].masks