# Google Drive raw data manifest — n1_repulsive

The generator raw file is `raw_macro_part000.h5` (118762744 bytes). Because the Drive connector did not accept the >100 MB local file reference, the exact HDF5 bytes were gzip-compressed for Drive transport only.

Drive transport file: `raw_macro_part000.h5.gz` (94178506 bytes)

- original HDF5 SHA-256: `b015efde42555f4b38046b9af0d64c35faf0f3f372f408f4fca5f1f025277465`
- gzip transport SHA-256: `e73c3c74fd4ae248e1a3b463c9f92b69cb68276b8ae7c50ac54b248942fdb5f7`

Reconstruction:

```sh
gunzip -c raw_macro_part000.h5.gz > raw_macro_part000.h5
```

The reconstructed HDF5 must match the original SHA-256 above. Compression changes no rows or state values.
