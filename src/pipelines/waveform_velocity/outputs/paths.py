"""Shared output-schema resolution for waveform-velocity serializers."""

from input_output.schema import EyeFlowOutputPaths


def resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)


__all__ = ["resolve_output_paths"]
