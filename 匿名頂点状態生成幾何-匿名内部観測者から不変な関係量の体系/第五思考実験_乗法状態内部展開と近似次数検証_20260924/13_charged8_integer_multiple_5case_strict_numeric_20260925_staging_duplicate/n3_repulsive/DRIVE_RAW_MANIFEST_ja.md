# Google Drive raw data manifest — n3_repulsive

The four generator HDF5 raw parts are preserved exactly. For Drive transport only, each HDF5 file was gzip-compressed and split into <=90 MiB chunks. Concatenating the chunks and gunzipping reconstructs the exact original HDF5 bytes.

## raw_macro_part000.h5

- original size: 325018498 bytes
- original SHA-256: `16d90d6973f4d5dfc307cd8f6d065ce14c3c2f82665ae0189437e1dd2133aa8b`
- gzip SHA-256: `a2e16bf1b40d9d44b8b6e0875e2d6f7d15869b1f7307e87f8cf4cec8a499e19e`
- Drive chunks:
  - `raw_macro_part000.h5.gz.part00` — 94371840 bytes — SHA-256 `20f314c3511f71ed0c0e1ffc19835cd63acd7c6efaef7800cc76985b9510d1c2`
  - `raw_macro_part000.h5.gz.part01` — 94371840 bytes — SHA-256 `107ee7ce18652c16518edc5cd80510355a72f42a9537b09311f915fe8377a34a`
  - `raw_macro_part000.h5.gz.part02` — 65424140 bytes — SHA-256 `0a5c56c72b88582f35acdd736a0ea909acd905b31b150e5fbd23b94295dee62a`

Reconstruction:

```sh
cat raw_macro_part000.h5.gz.part* > raw_macro_part000.h5.gz
gunzip -c raw_macro_part000.h5.gz > raw_macro_part000.h5
```

## raw_macro_part001.h5

- original size: 304275874 bytes
- original SHA-256: `69ad29ba76e7220aee01cc33d27ee4ff744398e9ba4e7c8e3ec297ac49954dea`
- gzip SHA-256: `59106a8e6d35ea2797d6e60b5d21f3e979e54b0a2ff0f8d4be40eb2785a3249d`
- Drive chunks:
  - `raw_macro_part001.h5.gz.part00` — 94371840 bytes — SHA-256 `e6954d375678d5d4141942d31ec74edbf6ce86cfab242c0edbef86cadb53ff51`
  - `raw_macro_part001.h5.gz.part01` — 94371840 bytes — SHA-256 `46e767de38bfa187961bb751035335d59b615e57bb65bd6ba9aaececf224f6c6`
  - `raw_macro_part001.h5.gz.part02` — 37684547 bytes — SHA-256 `d46b847c19efa00128bc605ab6d9579427fe8bd1f74386e889a6aab3524ce17c`

Reconstruction:

```sh
cat raw_macro_part001.h5.gz.part* > raw_macro_part001.h5.gz
gunzip -c raw_macro_part001.h5.gz > raw_macro_part001.h5
```

## raw_macro_part002.h5

- original size: 282688486 bytes
- original SHA-256: `481e148fb6da3dad0514d8d5a88dcf6f55dc06ce13962cbc1e378da3bff68811`
- gzip SHA-256: `24809b6dacdb4146ccf3b5d5dbd7b0628921fb1740894ed1dc52c21ff4746c64`
- Drive chunks:
  - `raw_macro_part002.h5.gz.part00` — 94371840 bytes — SHA-256 `38c43d0264854c62e458454e6c636a358aaca35acc52dafb3de34bb4d4b64e1f`
  - `raw_macro_part002.h5.gz.part01` — 94371840 bytes — SHA-256 `d01077ba9f4fba91d6828704838b43da53eb4467336b405044fb58c4436bba56`
  - `raw_macro_part002.h5.gz.part02` — 6277676 bytes — SHA-256 `8c8a4cef2165aa0a493bb9c3a13ef0a547c03e820179b0a5b45191f4b4fb31ad`

Reconstruction:

```sh
cat raw_macro_part002.h5.gz.part* > raw_macro_part002.h5.gz
gunzip -c raw_macro_part002.h5.gz > raw_macro_part002.h5
```

## raw_macro_part003.h5

- original size: 232481483 bytes
- original SHA-256: `13e0f4b6b443af94b65d0fbbf3f441b9aab5d3a6efd9758974e2c9ff3e776fc0`
- gzip SHA-256: `dd7c2ffe4d837e89605b3710c9e4c0d8c52fbf5dcb29a17b7745bc66ae4ea617`
- Drive chunks:
  - `raw_macro_part003.h5.gz.part00` — 94371840 bytes — SHA-256 `a964a7b9c1012430849a9b0963a83b71de6afc3112df31ac04ac6ae0efff5b23`
  - `raw_macro_part003.h5.gz.part01` — 69717198 bytes — SHA-256 `16e419a1665044410b61bc3755dba3f3d809b77980d57b0c17e2f4f7f80c96a6`

Reconstruction:

```sh
cat raw_macro_part003.h5.gz.part* > raw_macro_part003.h5.gz
gunzip -c raw_macro_part003.h5.gz > raw_macro_part003.h5
```

