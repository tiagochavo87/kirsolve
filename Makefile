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
	rm -rf _site
	mkdir -p _site/src/kirsolve
	cp docs/index.html _site/index.html
	cp app.py _site/app.py
	cp assets/logo.jpeg _site/logo.jpeg
	cp src/kirsolve/*.py _site/src/kirsolve/
	@echo "Agora: python -m http.server -d _site 8000"
	@echo "E abra http://localhost:8000"

app:
	streamlit run app.py

clean:
	rm -rf $(OUT) .pytest_cache **/__pycache__ dist build _site
