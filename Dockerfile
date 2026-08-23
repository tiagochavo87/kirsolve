FROM python:3.12-slim

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY pyproject.toml .
COPY src/ ./src/
COPY tests/ ./tests/
RUN pip install --no-cache-dir -e .

# dados de entrada e saida sao montados como volumes
VOLUME ["/data", "/out"]
ENTRYPOINT ["kirsolve"]
CMD ["--help"]
