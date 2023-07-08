from src.variants.ssim_multiscale_base import SSIMMultiScaleVariant
from src.variants.resized_ssim import ResizedSSIMVariant


class ResizedSSIMMultiScaleVariant(SSIMMultiScaleVariant, ResizedSSIMVariant):

    name = "Resized MultiScale Structural Similarity Index Measure"
