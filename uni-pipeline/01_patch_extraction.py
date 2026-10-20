
#!/usr/bin/env python3

import argparse
import os
import re

import cv2
import numpy as np
import pandas as pd


def sorted_key(file_path):
    """Extract a TMA identifier from a filepath."""
    match = re.search(r"(TMA\d+_\d+)", str(file_path))
    return match.group(1) if match else None


def remove_small_components(img, min_size=150):
    """Remove small connected components, or keep only the largest if min_size < 0."""
    nb_components, output, stats, _ = cv2.connectedComponentsWithStats(
        img.astype("uint8"), connectivity=8
    )
    sizes = stats[1:, -1]
    result = np.zeros(output.shape, dtype=np.uint8)

    if min_size < 0 and len(sizes):
        largest_component = np.argmax(sizes) + 1
        result[output == largest_component] = 255
    else:
        for i in range(nb_components - 1):
            if sizes[i] >= min_size:
                result[output == i + 1] = 255

    return result


def create_mask(dapi):
    """Create a binary tissue mask from an H&E image."""
    h, w = dapi.shape[:2]
    temp_img = dapi[::4, ::4]

    if temp_img.ndim == 3:
        temp_img = cv2.cvtColor(temp_img, cv2.COLOR_RGB2GRAY)

    temp_img = cv2.threshold(
        temp_img, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU
    )[1].astype("uint8")

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (15, 15))
    temp_img = cv2.morphologyEx(temp_img, cv2.MORPH_CLOSE, kernel)
    temp_img = remove_small_components(temp_img, min_size=-1)

    temp_img = cv2.resize(
        temp_img, (w, h), interpolation=cv2.INTER_NEAREST
    )
    return (temp_img > 128).astype("uint8")


def patching(image, mask, patch_size=224, stride=224):
    """Extract patches containing sufficient tissue."""
    patches, coords = [], []
    h, w = image.shape[:2]

    for y in range(0, h - patch_size + 1, stride):
        for x in range(0, w - patch_size + 1, stride):
            patch = image[y:y + patch_size, x:x + patch_size]
            mask_patch = mask[y:y + patch_size, x:x + patch_size]

            if mask_patch.mean() > 0.1:
                patches.append(patch)
                coords.append((x, y))

    return patches, coords


def patch_quality_control(
    patches, coords, min_tissue=0.1, min_color_std=2.0
):
    """Filter patches based on tissue content and color variation."""
    kept_patches, kept_coords = [], []

    for patch, coord in zip(patches, coords):
        if patch.mean() < 5:
            continue

        tissue = np.any(patch < 250, axis=-1)
        if tissue.mean() < min_tissue:
            continue

        color_std = patch.mean(axis=(0, 1)).std()
        if color_std < min_color_std:
            continue

        kept_patches.append(patch)
        kept_coords.append(coord)

    return kept_patches, kept_coords


def main(args):
    df = pd.read_csv(args.regions_csv)
    df = df.groupby("filename").first().reset_index()
    df = df[df["Tissue_region"].isin(["central", "meta", "peri"])]

    central_files = df["filename"].astype(str).tolist()
    paths_df = pd.read_csv(args.paths_csv)

    mapping = [
        path for path in paths_df["filepath"]
        if sorted_key(path) in central_files
    ]

    os.makedirs(args.output_dir, exist_ok=True)

    for path in mapping:
        image = np.load(path)
        mask = create_mask(image)

        masked_image = (
            image * mask[..., np.newaxis]
            + 255 * (1 - mask[..., np.newaxis])
        ).astype(image.dtype)

        patches, coords = patching(
            masked_image, mask, args.patch_size, args.stride
        )
        patches, coords = patch_quality_control(patches, coords)

        file_id = sorted_key(path)
        if not patches:
            continue

        np.save(
            os.path.join(args.output_dir, f"{file_id}_patches.npy"),
            np.array(patches),
        )
        np.save(
            os.path.join(args.output_dir, f"{file_id}_coords.npy"),
            np.array(coords),
        )

        print(f"{file_id}: saved {len(patches)} patches")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract and quality-control H&E tissue patches."
    )
    parser.add_argument("--regions-csv", required=True)
    parser.add_argument("--paths-csv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--patch-size", type=int, default=224)
    parser.add_argument("--stride", type=int, default=224)

    main(parser.parse_args())
