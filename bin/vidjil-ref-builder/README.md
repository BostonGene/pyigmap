## How to build vidjil reference

1. Build image
```bash
docker build -t vidjil-ref-builder .
```

2. Run
```bash
docker run --rm \
  -v ./:/tmp/ \ 
  vidjil-ref-builder \
  --output-archive /tmp/vidjil.germline.new.tar.gz
```
