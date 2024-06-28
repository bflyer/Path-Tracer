
def transform_string(input_str, x_mult=1, y_mult=1, z_mult=1, dx=0, dy=0, dz=0):
    # 将输入字符串分割成列表，每行为一个元素
    lines = input_str.strip().split('\n')
    
    transformed_lines = []
    for line in lines:
        parts = line.strip().split()
        if parts and parts[0] == 'v':  # 确保parts不为空且首元素为'v'
            # 提取坐标值并转换
            try:
                x, y, z = map(float, parts[1:])
                new_x = x * x_mult + dx
                new_y = y * y_mult + dy
                new_z = z * z_mult + dz
                transformed_lines.append(f"v {new_x} {new_y} {new_z}\n")
            except ValueError:
                print(f"警告: 无法解析坐标 '{line}', 跳过此行.")
        else:  # 对非 'v' 开头的行直接追加
            transformed_lines.append(line + '\n')  # 注意保持原有行结束符

    return "".join(transformed_lines)

def transform_file(input_filename, output_filename, x_mult=1, y_mult=1, z_mult=1, dx=0, dy=0, dz=0):
    try:
        with open(input_filename, 'r') as file:
            input_data = file.read()
    except FileNotFoundError:
        print(f"错误：找不到文件 {input_filename}")
        return
    
    transformed_data = transform_string(input_data, x_mult, y_mult, z_mult, dx, dy, dz)
    
    with open(output_filename, 'w') as file:
        file.write(transformed_data)
        
    print(f"处理完成，结果已保存至 {output_filename}")

# 定义输入输出文件名以及变换参数
input_filename = "input.txt"
output_filename = "output.txt"
x_multiplier = 50
y_multiplier = 50
z_multiplier = 50
x_translation = 50
y_translation = -6
z_translation = 80

# 调用函数处理文件
transform_file(input_filename, output_filename, x_multiplier, y_multiplier, z_multiplier, x_translation, y_translation, z_translation)