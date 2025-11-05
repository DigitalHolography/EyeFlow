from ultralytics import YOLO

def import_yolo_model(model_path):
    """
    Imports and returns a YOLO model from the specified path.

    Args:
        model_path: path to the YOLO segmentation model (e.g., YOLOv11)

    Returns:
        YOLO model instance
    """
    model = YOLO(model_path)
    return model