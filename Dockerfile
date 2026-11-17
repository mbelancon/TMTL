# Image to host the interface (Hugging Face Spaces, Render, a lab server…)
FROM python:3.12-slim
RUN pip install --no-cache-dir numpy scipy
WORKDIR /app
COPY tlfit ./tlfit
ENV TLFIT_HOST=0.0.0.0 PORT=7860 TLFIT_MAX_FITS=2 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
EXPOSE 7860
CMD ["python", "-m", "tlfit", "--no-browser"]
