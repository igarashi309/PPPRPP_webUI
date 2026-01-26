
import os
from PIL import Image, ImageDraw, ImageFont

def combine_figures(panel_paths, output_path, labels=["A", "B", "C"]):
    images = [Image.open(p) for p in panel_paths]
    
    # Calculate dimensions
    # We want uniform width, usually matches the widest
    max_width = max(img.width for img in images)
    
    # Resize images to match max_width (maintaining aspect ratio? No, usually browser/detail are wide enough)
    # Actually, if we resize, we might distort text. Better to paste onto a canvas centered or left aligned.
    # Let's assume user wants them left aligned or centered.
    
    total_height = sum(img.height for img in images)
    padding = 40
    header_height = 40
    
    # Create canvas
    canvas_height = total_height + (len(images) * (padding + header_height))
    canvas = Image.new('RGB', (max_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)
    
    # Try to load a font, or default
    try:
        font = ImageFont.truetype("Arial.ttf", 60)
    except:
        try:
             font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 60)
        except:
             font = ImageFont.load_default()

    current_y = 0
    for i, img in enumerate(images):
        # Draw Label (A, B, C)
        label = labels[i]
        draw.text((10, current_y), label, fill="black", font=font, stroke_width=2, stroke_fill="white")
        
        # Paste Image below label
        # Center the image?
        x_offset = (max_width - img.width) // 2
        canvas.paste(img, (x_offset, current_y + header_height))
        
        current_y += img.height + header_height + padding
        
    canvas.save(output_path)
    print(f"Created {output_path}")

# Figure 1
base_dir = "/Users/sola309/Downloads/PPPRPP_webUI-main/manuscript/figures"
combine_figures(
    [os.path.join(base_dir, f) for f in ["Fig1A_Motif.png", "Fig1B_Browser.png", "Fig1C_Detail.png"]],
    os.path.join(base_dir, "Figure_1.png")
)

# Figure S1
combine_figures(
    [os.path.join(base_dir, f) for f in ["FigS1A_Motif.png", "FigS1B_Browser.png", "FigS1C_Detail.png"]],
    os.path.join(base_dir, "Figure_S1.png")
)
