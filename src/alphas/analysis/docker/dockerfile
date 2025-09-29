# Use Debian-slim image with python 3.10
# Note: Use Bookworm to get gcc12 since gcc14 barfs when compiling gvar
# Note: Use python 3.10 required for statistics package

#---------- builder

FROM python:3.10-slim-bookworm AS builder

# Install gcc
RUN apt update \
    && \
    apt install -y \
    build-essential \
    libhdf5-dev \
    && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory inside the container
WORKDIR /app

# Copy only the requirements file first to leverage Docker's layer caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --upgrade pip && pip install --no-cache-dir -r requirements.txt

#---------- final

FROM python:3.10-slim-bookworm AS final

RUN apt update \
    && \
    apt install -y \
    wget \
    libhdf5-103 \
    && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory inside the container
WORKDIR /app

# Copy only the requirements file first to leverage Docker's layer caching
COPY requirements.txt .

# Copy site  python packages
COPY --from=builder /usr/local/lib/python3.10/site-packages /usr/local/lib/python3.10/site-packages

RUN pip install --upgrade pip && pip install --no-cache-dir streamlit

# Copy the rest of the application code
COPY data_app.py plot.py ensemble_not_found.png .

# Expose the port your web application listens on
EXPOSE 8501

# Define the command to run your application
CMD ["streamlit", "run",  "--server.address", "0.0.0.0", "--server.port=8501", "data_app.py"]

