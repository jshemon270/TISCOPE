TISCOPE: TISsue Spatial Characterization and Omics PipelinE

Screenshot 2026-10-19 at 16 26 35
01. Cellpose-Based Cell Segmentation
Cell segmentation is performed using Cellpose (Stringer et al., 2021) on the DAPI channel. In our use case, cyclic immunofluorescence images are stored in a NumPy array, with the first DAPI image accessed using image[0, 0, :, :].

The resulting segmentation masks and the corresponding intensity images are saved as separate NumPy arrays.

Segmentation masks

The DAPI-derived cell segmentation masks are saved as TMA_masks.npy:

np.save(os.path.join(save_dir, f"{base}_masks.npy"), masks)
Intensity images

The intensity images for the remaining staining channels (excluding DAPI) are saved as TMA_intensity.npy:

np.save(os.path.join(save_dir, f"{base}_intensity.npy"), intensity_images)
These outputs are used in downstream analyses, including cellular feature extraction and spatial characterization.

02. Segmentation Mask Extraction
This notebook saves several masks derived from the DAPI mask and the original intensity images using scikit-image (van der Walt et al., 2014). These are all numpy array d.types:

1.) An overlap mask, overlapping the DAPI (most cells, least tissue degradation) to the final cycle (least cells, most tissue degradation)

overlap = (DAPI_file & (b_last_cycle_DAPI > 0)).astype(bool)
It then performs a watershed on the overlapped image.

2.) Following the watershed, the nuclei are expanded by 15 pixels, assuming our cell body mask.

expanded_overlap = expand_labels(seg1, distance=15)
3.) We then save another mask that holds just the cell body, without the nuclei.

cell_body_mask = (expanded_overlap > 0) & (overlap == 0)
The final two masks are the DAPI mask (again) and the intensity images (again), so they stay in the proper folder.

np.save(os.path.join(output_dir, f"{file_prefix}_cell_body_mask.npy"), cell_body_mask)
np.save(os.path.join(output_dir, f"{file_prefix}_expanded_overlap.npy"), expanded_overlap)
np.save(os.path.join(output_dir, f"{file_prefix}_overlap_mask.npy"), overlap)
np.save(os.path.join(output_dir, f"{file_prefix}_dapi_masks.npy"), DAPI_file)
# Slicing the intensity image to exclude DAPI cycle and save space
np.save(os.path.join(output_dir, f"{file_prefix}_int_img_wo_DAPI.npy"), Int_Img_File[:, 1:, :, :])
03. Region Property Extraction
The region prop extraction uses the overlap image as the cell label and the intensity image as either the raw DAPI or the other cycles images. This function is also contained in scikit-image.

props = regionprops_table(
        labels,
        intensity_image= raw_dapi,
        properties=[
            'label', 'centroid', 'orientation',
            'axis_major_length', 'axis_minor_length',
            'area', 'eccentricity', 'perimeter', 'mean_intensity'
        ]
    )
This extracts all relevant features for further analysis.

04. UNIform Image Normalization
We chose UNIform (Wang et al., 2025) normalization for our mean intensities as its specifically designed for multiplex images. Here is an example output below: image

Derived from here this paper: Link Text

05. K-Means Clustering and UMAP Visualization
Kmeans clustering (Pedregosa et al., 2011) was done to delegate clusters in our dataset, a silhouette score was used to assess. A UMAP (McInnes et al., 2018) was then performed to analyze it (cuml.UMAP improves speed substantially, not used here).

06. Z-Score-Based Feature Analysis
A Z-score was used to see what features deviate from the mean within these clusters. Establishing our relevant biomarkers within each cluster. Using this we labelled the clusters accordingly.

07. Spatial Metrics
These metrics assess relevant spatial analysis between these clusters. We evaluated spatial relationships using neighborhood enrichment, centrality, co_occurrence, and ripley scores. All derived from squidpy (Palla et al., 2022): Link Text

08. Spatial Clustering Analysis
Following CellCharter (Varrone et al., 2024) documentation I applied spatial clusters to the dataset as well, incorporating the positional information, resulting in 3 spatial domains. A Gaussian Mixture Model was used to determine the K, and BIC scores were used the assess. They were validated using stability scores and visualized with a elbow curve.

The link to the documentation: Link Text

We then used the Z-score from above to test what biomarkers are present in the clusters, this led to 3 broader categories.

09. Cliff's Delta Effect Size Analysis
Cliffs delta measures the magnitude of difference between two groups, it can be used in tandem with statistical tests. While statistical tests test against a null hypothesis, Cliffs delta simply tells how strongly two groups differ.

10. XGBoost Classification
XGBoost (Chen & Guestrin, 2016) was used to predict the TNM staging system.

11. Pair Correlation Function (PCF) Analysis
12. Whole-Tissue Pair Correlation Function Analysis
13. Clark–Evans Spatial Clustering Analysis
