## How to build archive with OLGA models

1. Build image
```bash
docker build -t olga-ref-builder .
```

2. Run
```bash
docker run --rm \
  -v ./:/tmp/ \ 
  olga-ref-builder \
  --out-archive /tmp/olga-models.tar.gz
```