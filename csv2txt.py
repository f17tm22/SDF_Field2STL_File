import csv
 
# 读取CSV文件
csv_file_path = 'dpp.csv'
txt_file_path = 'dpp_txt.txt'
 
with open(csv_file_path, 'r') as csv_file, open(txt_file_path, 'w') as txt_file:
    # 创建CSV读取器
    csv_reader = csv.reader(csv_file)
 
    # 逐行读取CSV文件，将每行的内容以制表符分隔写入txt文件
    for row in csv_reader:
        txt_file.write(' '.join(row) + '\n')
 
print(f"Successfully converted {csv_file_path} to {txt_file_path} using csv module.")
