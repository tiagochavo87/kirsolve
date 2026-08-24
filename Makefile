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

site:
	pip install build
	python -m build --wheel
	cp dist/kirsolve-*-py3-none-any.whl docs/
	cp app.py assets/logo.jpeg docs/
	@echo "Agora: python -m http.server -d docs 8000"
	@echo "E abra http://localhost:8000"

app:
	streamlit run app.py

clean:
	rm -rf $(OUT) .pytest_cache **/__pycache__ dist build
	rm -f docs/*.whl docs/app.py docs/logo.jpeg
