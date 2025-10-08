FROM rootproject/root:6.28.04-ubuntu22.04

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .
RUN make

CMD ["./likelihood_fit"]
