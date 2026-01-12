document.addEventListener('mousemove', function (e) {
    console.log('鼠标位置：', e.clientX, e.clientY);
    // 获取鼠标当前位置
    const mouseX = e.clientX;
    const mouseY = e.clientY;

    // 创建一个小圆点元素
    const circle = document.createElement('div');
    circle.classList.add('circle');

    // 随机生成一个方向（例如：X轴和Y轴方向的偏移量）
    const randomAngle = Math.random() * 360;  // 0 到 360 度之间的随机角度
    const randomDistance = Math.random() * 50;  // 随机距离，控制圆点偏移量
    const randomX = randomDistance * Math.cos(randomAngle);
    const randomY = randomDistance * Math.sin(randomAngle);

    // 设置圆点的初始位置
    circle.style.left = `${mouseX + randomX}px`;
    circle.style.top = `${mouseY + randomY}px`;

    // 将小圆点添加到页面
    document.getElementById('container').appendChild(circle);

    // 在1.5秒后移除小圆点
    setTimeout(() => {
        circle.remove();
    }, 1500);  // 1500ms即1.5秒后删除
});