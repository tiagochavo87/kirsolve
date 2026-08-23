IMAGE ?= kirsolve:0.1.0
INPUT ?= data/entrada.xlsx
OUT   ?= resultados

.PHONY: install test run docker docker-run clean

install:
	pip install -e .

test:
	pytest -q

run:
	python -m kirsolve run $(INPUT) -o $(OUT) --sheet-tool Sheet5=kir-mapper

inspect:
	python -m kirsolve inspect $(INPUT)

docker:
	docker build -t $(IMAGE) .

docker-run:
	docker run --rm -v $(PWD):/work -w /work $(IMAGE) \
		run $(INPUT) -o $(OUT) --sheet-tool Sheet5=kir-mapper

clean:
	rm -rf $(OUT) .pytest_cache **/__pycache__
