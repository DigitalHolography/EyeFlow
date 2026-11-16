"""Resolve the run's authoritative annulus geometry in native image pixels."""


def resolve_annulus_geometry(ctx, optic_disc, image_shape):
    """Use intentional EyeFlow settings, never legacy acquisition attributes.

    Only a count explicitly loaded from a saved EyeFlow configuration requests
    legacy non-overlapping sampling. All consumers reuse this prepared geometry.
    """
    settings = getattr(ctx, "annulus_settings", {})
    return optic_disc.annulus_geometry(
        image_shape,
        annulus_size=settings.get("annulus_size"),
        annulus_stride=settings.get("annulus_stride"),
        number_of_radii_in_fov=settings.get("number_of_radii_in_FOV"),
    )
