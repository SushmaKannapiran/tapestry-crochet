import sys
sys.path.insert(0, '/Users/ksushma/Desktop/TB Cursor/pylibs')
from fpdf import FPDF

class PDF(FPDF):
    def header(self):
        pass

    def section_title(self, title):
        self.set_font('Helvetica', 'B', 16)
        self.set_text_color(30, 64, 175)
        self.cell(0, 10, title, new_x="LMARGIN", new_y="NEXT")
        self.set_draw_color(30, 64, 175)
        self.line(10, self.get_y(), 200, self.get_y())
        self.ln(4)

    def sub_title(self, title):
        self.set_font('Helvetica', 'B', 12)
        self.set_text_color(55, 65, 81)
        self.cell(0, 8, title, new_x="LMARGIN", new_y="NEXT")
        self.ln(2)

    def body_text(self, text, bold=False):
        self.set_font('Helvetica', 'B' if bold else '', 10)
        self.set_text_color(26, 26, 26)
        self.multi_cell(0, 6, text)
        self.ln(2)

    def add_table(self, headers, data, col_widths=None):
        if col_widths is None:
            col_widths = [190 / len(headers)] * len(headers)

        self.set_font('Helvetica', 'B', 9)
        self.set_fill_color(30, 64, 175)
        self.set_text_color(255, 255, 255)
        for i, h in enumerate(headers):
            self.cell(col_widths[i], 8, h, border=1, fill=True, align='C')
        self.ln()

        self.set_font('Helvetica', '', 8.5)
        self.set_text_color(26, 26, 26)
        fill = False
        for row in data:
            if fill:
                self.set_fill_color(248, 250, 252)
            else:
                self.set_fill_color(255, 255, 255)
            max_h = 8
            for i, cell in enumerate(row):
                self.cell(col_widths[i], max_h, str(cell), border=1, fill=True, align='L')
            self.ln()
            fill = not fill
        self.ln(4)

    def bullet(self, text):
        self.set_font('Helvetica', '', 10)
        self.set_text_color(26, 26, 26)
        self.cell(6, 6, '-')
        self.multi_cell(0, 6, text)
        self.ln(1)

    def numbered_item(self, num, text):
        self.set_font('Helvetica', 'B', 10)
        self.set_text_color(30, 64, 175)
        self.cell(8, 6, f"{num}.")
        self.set_font('Helvetica', '', 10)
        self.set_text_color(26, 26, 26)
        self.multi_cell(0, 6, text)
        self.ln(1)


pdf = PDF()
pdf.set_auto_page_break(auto=True, margin=15)
pdf.add_page()

# Title
pdf.set_font('Helvetica', 'B', 24)
pdf.set_text_color(26, 26, 26)
pdf.cell(0, 14, 'India Crochet & Yarn Market Research', new_x="LMARGIN", new_y="NEXT")
pdf.set_font('Helvetica', '', 14)
pdf.set_text_color(107, 114, 128)
pdf.cell(0, 10, 'The REAL Numbers - Hands Dirty Edition', new_x="LMARGIN", new_y="NEXT")
pdf.set_font('Helvetica', '', 11)
pdf.cell(0, 8, 'Date: March 2026', new_x="LMARGIN", new_y="NEXT")
pdf.ln(6)
pdf.set_draw_color(229, 231, 235)
pdf.line(10, pdf.get_y(), 200, pdf.get_y())
pdf.ln(8)

# Section 1
pdf.section_title('1. YouTube: Indian Knitting/Crochet Channels')
pdf.body_text('Source: Modash Influencer Database (Last updated: March 14, 2026)')
pdf.body_text('Total found: 44 knitting/crochet YouTubers in India')
pdf.sub_title('Top 14 Channels by Subscribers')

yt_headers = ['#', 'Channel', 'Subscribers', 'Top Video Views', '% India']
yt_data = [
    ['1', 'Ritu Creations', '883K', '2.5M views', '75%'],
    ['2', 'Magical Threadz (Rajnish)', '730K', '4.8M views', '82%'],
    ['3', "Chinkiz Knitting Knife", '586K', '11.8M views', '57%'],
    ['4', 'Creativity Lovers', '487K', '2M views', '88%'],
    ['5', 'PIHOO TV', '396K', '4.8M views', '70%'],
    ['6', 'Housie Things', '298K', '28.3M views', 'N/A'],
    ['7', 'Anita Dekhe Aur Sikhe', '183K', '2.4M views', '79%'],
    ['8', 'Bandna Bhardwaj Vlogs', '136K', '1.2M views', '93%'],
    ['9', 'Homemade Crafting (Arti)', '126K', '948K views', '54%'],
    ['10', "Krishnas Creative Knitting", '115K', '2.5M views', '69%'],
    ['11', 'Sanwariya Crafts', '56.9K', '825K views', '97%'],
    ['12', 'Pooja Creations', '46.3K', '717K views', '97%'],
    ['13', 'Knit With Ammu (Malayalam)', '17.5K', '109K views', '66%'],
    ['14', 'INDI DIY (North-East)', '10.8K', '161K views', '37%'],
]
pdf.add_table(yt_headers, yt_data, [10, 62, 30, 38, 50])
pdf.body_text('Total subscribers across top 14 channels: ~4.07 million', bold=True)
pdf.body_text('Key Insight: This is just the top 14 out of 44 found. Many channels are Hindi-language traditional knitting (sweaters, Laddu Gopal dresses). The modern crochet/amigurumi wave is a separate, newer segment.')

# Section 2
pdf.section_title('2. Instagram: Indian Crochet Accounts')
ig_headers = ['Account', 'Followers', 'Type']
ig_data = [
    ['@magicneedlesindia', '126K', 'Yarn store + brand'],
    ['@crochet_by_kunal', '124K', 'Creator + shop'],
    ['Crochet Now India', 'Est. 20-30K', 'Yarn store (Hyderabad)'],
    ['Various small creators', '1K-50K each', 'Hundreds of them'],
]
pdf.add_table(ig_headers, ig_data, [60, 50, 80])
pdf.body_text('Global benchmark: The #crochet hashtag has 26.7 million posts on Instagram globally. India-specific hashtags (#crochetindia, #indiancrocheters) are a growing fraction.')

# Section 3
pdf.section_title('3. Reddit')
pdf.body_text('r/IndianHobbyists - searched for "crochet" - ZERO results.')
pdf.body_text('The crochet community in India is NOT on Reddit. They are on:')
pdf.bullet('Facebook groups - biggest platform for Indian crocheters')
pdf.bullet('Instagram - creators + sellers')
pdf.bullet('YouTube - tutorials (Hindi, Malayalam, English)')
pdf.bullet('WhatsApp groups - e.g., MICQ operates via WhatsApp')

# Section 4
pdf.section_title('4. Online Yarn Stores - Customer Data')
store_headers = ['Store', 'Claimed Customers', 'Based in', 'Est. Since']
store_data = [
    ['Magic Needles', '300,000+ lifetime', 'Mumbai', '2015'],
    ['Crochet Now India', '20,000+', 'Hyderabad', '2017'],
    ['Pony Craft Store', 'Unknown', 'Unknown', '2012'],
    ['ABC Wools', 'Unknown', 'Rajasthan', 'Unknown'],
    ['WoolCraft Store', 'Unknown', 'Unknown', 'Unknown'],
    ['Knitting Happiness', 'Unknown', 'Unknown', 'Since 1954'],
    ['Avya Store', 'Unknown', 'Unknown', 'Unknown'],
    ['Shuttles & Needles', 'Unknown', 'Unknown', 'Unknown'],
]
pdf.add_table(store_headers, store_data, [50, 45, 50, 45])
pdf.body_text('Note: Amazon India and Flipkart also sell significant volumes of yarn - uncounted buyers.')

# Section 5
pdf.section_title('5. Viral Moments - Demand Proof')
viral_headers = ['What', 'Numbers', 'When']
viral_data = [
    ["Nikita Shaw sunflower clutcher reel", '33 million views', '2025'],
    ['Stushe crochet bags D2C brand', 'Rs 1 Cr in < 1 year', '2024-25'],
    ["Chinkiz Knitting Knife top video", '11.8 million views', '-'],
    ['Housie Things top video', '28.3 million views', '-'],
    ['Magical Threadz crochet jacket', '4.8 million views', '-'],
    ['MICQ Guinness World Records', '6 records', '2017-2025'],
]
pdf.add_table(viral_headers, viral_data, [85, 55, 50])

# Section 6
pdf.add_page()
pdf.section_title('6. Market Size Estimate - Bottom-Up')
pdf.sub_title('Method: Bottom-up from actual touchpoints')
est_headers = ['Signal', 'Raw Number', 'Logic', 'Est. People']
est_data = [
    ['YouTube subs (top 14)', '4.07M subs', '~30% actually craft', '~1.2M'],
    ['YouTube total (44+ ch)', '5-6M total subs', 'Dedup, 30% active', '~1.5M'],
    ['Magic Needles customers', '300K lifetime', 'x5-8 (many other stores)', '1.5-2.4M'],
    ['Traditional knitters', 'Large - no online', 'Kashmir, Himachal, UP', '2-5M'],
    ['Instagram followers', '300-400K combined', '50% overlap w/ YT', '~200K add.'],
    ['Facebook groups', '6K+ MICQ + locals', 'Many small groups', '50-100K add.'],
    ['Offline only', 'Unknown', 'Small towns, traditional', '1-3M'],
]
pdf.add_table(est_headers, est_data, [42, 38, 55, 55])

# Section 7
pdf.section_title('7. The Bottom Line')
bl_headers = ['Segment', 'Estimate', 'Confidence']
bl_data = [
    ['Active hobbyist crocheters', '300K - 500K', 'Medium-High'],
    ['Active knitters (traditional)', '1.5 - 3M', 'Medium'],
    ['Total active yarn crafters', '2 - 3.5M', 'Medium'],
    ['Casual / occasional', '1 - 2M', 'Low-Medium'],
    ['Following/interested', '3 - 5M', 'Low'],
]
pdf.add_table(bl_headers, bl_data, [80, 55, 55])

pdf.set_font('Helvetica', 'B', 12)
pdf.set_text_color(30, 64, 175)
pdf.cell(0, 10, 'Key Numbers:', new_x="LMARGIN", new_y="NEXT")
pdf.set_font('Helvetica', 'B', 11)
pdf.set_text_color(26, 26, 26)
pdf.bullet('Active crocheters in India: ~300K - 500K')
pdf.bullet('Active yarn crafters total (knit + crochet): ~2 - 3.5 million')

# Section 8
pdf.ln(4)
pdf.section_title('8. Why It Still Matters')
pdf.numbered_item(1, 'Market value: 300K-500K crocheters x Rs 2K-5K/year = Rs 60-250 crore market for crochet yarn alone.')
pdf.numbered_item(2, 'Growth rate: 20-30% year-on-year (post-COVID explosion continues).')
pdf.numbered_item(3, 'Conversion: Traditional knitting base (2-3M) can be converted to crochet.')
pdf.numbered_item(4, 'Proof: Magic Needles built Rs 8.5 Cr from Rs 1 lakh - just ONE brand.')
pdf.numbered_item(5, 'Global: Export multiplies TAM by 10x+. US has 45M knitters/crocheters.')

# Section 9
pdf.section_title('9. Competitive Landscape - Indian Yarn Brands')
comp_headers = ['Brand', 'Colors', 'Price/100g', 'Mercer.', 'OEKO']
comp_data = [
    ['Hobby Store (Magic Needles)', '61', 'Rs 155', 'Yes', 'Yes'],
    ['Vardhman', '~31', 'Rs 80-160', 'No', 'No'],
    ['Ayush Crafts (Urshii)', '~28', 'Rs 150-200', 'Yes', 'No'],
    ['Ganga', '~25', 'Rs 110-220', 'No', 'No'],
    ['Loombastic Studios', '~14', 'Rs 150+', 'Yes', 'No'],
    ['Hilush (Bangalore)', '7', 'Rs 200+', 'No', 'No'],
]
pdf.add_table(comp_headers, comp_data, [58, 22, 40, 35, 35])

pdf.sub_title('International Benchmark')
intl_headers = ['Brand', 'Colors', 'Country', 'Price/50g', 'Made in']
intl_data = [
    ['Scheepjes Catona', '113', 'Netherlands', 'Rs 165-250', 'INDIA'],
    ['Ricorumi (Rico)', '80', 'Germany', 'Rs 200+', 'China'],
    ['Paintbox Cotton DK', '~80', 'UK', 'Rs 200+', 'Unknown'],
    ['Lion Brand 24/7', '54', 'USA', 'Rs 500+', 'Unknown'],
]
pdf.add_table(intl_headers, intl_data, [42, 22, 36, 42, 48])
pdf.body_text("The gap: Indias best brand has 61 colors. Worlds best has 113 - manufactured in India. Nobody has 150+.", bold=True)

# Section 10
pdf.section_title('10. The Opportunity')
opp_headers = ['Factor', 'Detail']
opp_data = [
    ['Manufacturing', 'GP Yarns (Tirupur) offers 15,000 shades'],
    ['Raw material', "India is worlds #1 cotton producer"],
    ['Production cost', 'Rs 20-45 per 50g ball (mercerised cotton)'],
    ['Retail price', 'Rs 120-160/50g (India) or $1.50-2/50g (export)'],
    ['Color range target', '150+ colors = #1 in world for cotton craft yarn'],
    ['Domestic market', 'Rs 60-250 crore, growing 20-30% YoY'],
    ['Global market', '$4.2B crochet market, 45M crafters in US'],
    ['Startup benchmark', 'Magic Needles: Rs 1 lakh to Rs 8.5 Cr revenue'],
]
pdf.add_table(opp_headers, opp_data, [45, 145])

# Footer
pdf.ln(8)
pdf.set_draw_color(229, 231, 235)
pdf.line(10, pdf.get_y(), 200, pdf.get_y())
pdf.ln(4)
pdf.set_font('Helvetica', 'I', 8)
pdf.set_text_color(107, 114, 128)
pdf.multi_cell(0, 5, 'Research compiled March 2026. Sources: Modash, Scheepjes.com, MagicNeedles.in, CrochetNow.in, Yarnalicious, GP Yarns Tirupur, MoneyMint, Guinness World Records, Craft Yarn Council, various Instagram/YouTube accounts.')

pdf_path = '/Users/ksushma/Desktop/India_Crochet_Market_Research.pdf'
pdf.output(pdf_path)
print(f'PDF saved to: {pdf_path}')
