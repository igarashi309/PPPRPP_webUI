
import os
from PIL import Image, ImageDraw, ImageFont

def create_pdf_figure(panel_paths, output_path, labels=["A", "B", "C"]):
    # A4 Size at 300 DPI
    A4_WIDTH = 2480
    A4_HEIGHT = 3508
    DPI = 300
    
    # Create Canvas
    canvas = Image.new('RGB', (A4_WIDTH, A4_HEIGHT), "white")
    draw = ImageDraw.Draw(canvas)
    
    # Font Setup
    font_size = 100
    font_files = [
        "/Library/Fonts/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/System/Library/Fonts/Supplemental/Arial.ttf"
    ]
    font = None
    for f in font_files:
        if os.path.exists(f):
            try:
                font = ImageFont.truetype(f, font_size)
                print(f"Using font: {f}")
                break
            except Exception as e:
                print(f"Failed to load {f}: {e}")
                
    if font is None:
        font = ImageFont.load_default()
        # Scale for default font is harder, hopefully we found a ttf
        print("Using default font (might be small)")

    # Layout Parameters
    margin_left = 150
    margin_top = 200
    spacing = 100 # Space between panel and previous image
    label_offset = 20 # Space between Label and Image top
    
    current_y = margin_top
    
    for i, p_path in enumerate(panel_paths):
        if not os.path.exists(p_path):
            print(f"Warning: {p_path} not found.")
            continue
            
        img = Image.open(p_path)
        
        # Draw Label
        label = labels[i]
        
        # Calculate text bounding box to align
        # draw.textbbox is available in newer Pillow, else manual
        draw.text((margin_left, current_y), label, fill="black", font=font)
        
        # Image Position
        # We assume image label is top-left of the image block.
        # Actually standard is Label at Top-Left, Image below or beside.
        # Let's put Image below the label with slight offset.
        
        img_y = current_y + font_size + label_offset
        
        # Centering image or Left Align?
        # If we align left with margin:
        img_x = margin_left
        
        # Check if image needs scaling to fit?
        # A4 width 2480. Margin 150*2 = 300. Max width ~2180.
        # Image width 1470. Fits easily.
        # Optional: Upscale to fill width? 
        # Better to keep 1:1 to avoid blur, or modest upscale. 1470 is good size.
        
        canvas.paste(img, (img_x, img_y))
        
        # Update Y
        current_y = img_y + img.height + spacing

    # Save as PDF
    # PDF specific save
    canvas.save(output_path, "PDF", resolution=DPI)
    print(f"Saved {output_path}")

base_dir = "/Users/sola309/Downloads/PPPRPP_webUI-main/manuscript/figures"

# Figure 1
create_pdf_figure(
    [os.path.join(base_dir, f) for f in ["Fig1A_Motif.png", "Fig1B_Browser.png", "Fig1C_Detail.png"]],
    os.path.join(base_dir, "Figure_1.pdf")
)

# Figure S1
create_pdf_figure(
    [os.path.join(base_dir, f) for f in ["FigS1A_Motif.png", "FigS1B_Browser.png", "FigS1C_Detail.png"]],
    os.path.join(base_dir, "Figure_S1.pdf")
)
